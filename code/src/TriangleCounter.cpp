//
// Created by Cristian Boldrin on 17/11/25.
//

#include <TriangleCounter.h>

TriangleCounter::TriangleCounter(size_t aux_head_sample_size, size_t aux_tail_sample_size,
                                 size_t head_budget, size_t tail_budget, size_t random_seed) :
        sample_head_size_(head_budget),
        sample_tail_size_(tail_budget),
        aux_sample_size_head_(aux_head_sample_size),
        aux_sample_size_tail_(aux_tail_sample_size),
        rand_gen_(random_seed),
        aux_head_rand_indexer_(aux_head_sample_size, rand_gen_.getEngine()),
        aux_tail_rand_indexer_(aux_tail_sample_size, rand_gen_.getEngine()),
        head_rand_indexer_(head_budget, rand_gen_.getEngine()),
        tail_rand_indexer_(tail_budget, rand_gen_.getEngine()) {

    printf("\n----- Initializing TriangleCounter... -----\n");

    // -- init samples
    // node head and tail sample
    head_nodes_ = my_hash_set_t();
    tail_nodes_ = my_hash_set_t();
    // auxiliary edge sample
    aux_sample_head_ = std::vector<edge_t>(aux_sample_size_head_);
    aux_sample_tail_ = std::vector<edge_t>(aux_sample_size_tail_);
    // edge head and tail sample
    head_sample_ = std::vector<edge_t>(sample_head_size_);
    tail_sample_ = std::vector<edge_t>(sample_tail_size_);

    // -- reserve budgets
    subgraph_head_.reserve(2 * sample_head_size_);
    aux_subgraph_head_.reserve(2 * aux_sample_size_head_);
    subgraph_tail_.reserve(2 * sample_tail_size_);
    aux_subgraph_tail_.reserve(2 * aux_sample_size_tail_);

    printf("Initialized |A_H| =  %zu, |A_T| = %zu - |M_H| = %zu - |M_T| = %zu edges.\n",
           aux_sample_size_head_, aux_sample_size_tail_, sample_head_size_, sample_tail_size_);


}

TriangleCounter::~TriangleCounter() {
    // -- clear samples
    aux_sample_head_.clear();
    aux_sample_tail_.clear();
    head_sample_.clear();
    tail_sample_.clear();
    subgraph_head_.clear();
    subgraph_tail_.clear();
    aux_subgraph_head_.clear();
    aux_subgraph_tail_.clear();
    local_triangle_head_cnt_.clear();
    local_triangle_tail_cnt_.clear();
}


/**
 * Add edge to underlying subgraph (edge sample) - head, tail or aux
 * @param u
 * @param v
 * @param graph
 */
void TriangleCounter::add_edge(int u, int v, graph_t &graph) {
    // -- add edge to graph
    auto it_u = graph.find(u);
    if (it_u != graph.end()) it_u->second.emplace(v);
    else {
        emhash8::HashSet<int> neighs;
        neighs.emplace(v);
        graph.emplace_unique(u, std::move(neighs));
    }

    auto it_v = graph.find(v);
    if (it_v != graph.end()) {
        it_v->second.emplace(u);
    } else {
        emhash8::HashSet<int> neighs;
        neighs.emplace(u);
        graph.emplace_unique(v, std::move(neighs));
    }
}

/**
 * Delete edge from underlying subgraph (edge sample) - head, tail or aux
 * @param u
 * @param v
 * @param graph
 * @return
 */
bool TriangleCounter::delete_edge(int u, int v, graph_t &graph) {

    // -- delete v from u's neighbor list
    auto it_u = graph.find(u);
    if (it_u == graph.end()) return false;
    if (it_u->second.erase(v) == 0) return false;
    if (it_u->second.empty()) graph.erase(it_u);

    // -- delete u from v's neighbor list
    auto it_v = graph.find(v);
    if (it_v == graph.end()) return false;
    if (it_v->second.erase(u) == 0) return false;
    if (it_v->second.empty()) graph.erase(it_v);

    return true;

}


/**
 * Count triangles formed by edge {u, v} in the given subgraph and auxiliary subgraph
 * @param u
 * @param v
 * @param set
 * @param subgraph
 * @param aux_subgraph
 * @param cur
 * @param size
 * @param aux_cur
 * @param aux_size
 * @param local_triangles
 * @return
 */
double TriangleCounter::count_triangles(int u, int v, bool is_head, const my_hash_set_t &set,
                                        const graph_t &subgraph, const graph_t &aux_subgraph,
                                        size_t cur, size_t size, size_t aux_cur, size_t aux_size,
                                        node_triangle_t &local_triangles) {

    const emhash8::HashSet<int> *u_neighs = nullptr, *v_neighs = nullptr;
    const emhash8::HashSet<int> *u_neighs_aux = nullptr, *v_neighs_aux = nullptr;

    if (auto it = subgraph.find(u); it != subgraph.end()) u_neighs = &it->second;
    if (auto it = aux_subgraph.find(u); it != aux_subgraph.end()) u_neighs_aux = &it->second;
    if (auto it = subgraph.find(v); it != subgraph.end()) v_neighs = &it->second;
    if (auto it = aux_subgraph.find(v); it != aux_subgraph.end()) v_neighs_aux = &it->second;

    // Compute degrees
    int du = 0, dv = 0;
    if (u_neighs) du += u_neighs->size();
    if (u_neighs_aux) du += u_neighs_aux->size();
    if (v_neighs) dv += v_neighs->size();
    if (v_neighs_aux) dv += v_neighs_aux->size();

    // Swap to iterate through smaller degree vertex
    if (du > dv) {
        std::swap(u, v);
        std::swap(u_neighs, v_neighs);
        std::swap(u_neighs_aux, v_neighs_aux);
    }

    // Cache set membership for u and v (these are checked for EVERY triangle)
    const bool u_in_set = set.find(u) != set.end();
    const bool v_in_set = set.find(v) != set.end();

    // Precompute increments (guard underflow) using correction factors from reservoir sampling
    double inc_mm = 1.0;
    if (cur >= 2 && size >= 2) {
        inc_mm = (double(cur) * double(cur - 1)) / (double(size) * double(size - 1));
        if (inc_mm < 1.0) inc_mm = 1.0;
    }

    double inc_aa = 1.0;
    if (aux_cur >= 2 && aux_size >= 2) {
        inc_aa = (double(aux_cur) * double(aux_cur - 1)) / (double(aux_size) * double(aux_size - 1));
        if (inc_aa < 1.0) inc_aa = 1.0;
    }

    double inc_ma = 1.0;
    if (cur >= 1 && aux_cur >= 1 && size >= 1 && aux_size >= 1) {
        inc_ma = (double(cur) / double(size)) * (double(aux_cur) / double(aux_size));
        if (inc_ma < 1.0) inc_ma = 1.0;
    }

    double cum_count = 0.0;

    // Helper macro to update triangle counts efficiently
#define UPDATE_TRIANGLES(w, increment) do { \
        const bool w_in_set = set.find(w) != set.end(); \
        const int num_nodes = (int)(u_in_set + v_in_set + w_in_set); \
        const double contrib = (increment) * num_nodes; \
        cum_count += contrib; \
        \
        if (!(is_head && !u_in_set)) { \
            auto it_u = local_triangles.find(u); \
            if (it_u != local_triangles.end()) it_u->second += (increment); \
            else local_triangles.emplace_unique(u, (increment)); \
        } \
        \
        if (!(is_head && !v_in_set)) { \
            auto it_v = local_triangles.find(v); \
            if (it_v != local_triangles.end()) it_v->second += (increment); \
            else local_triangles.emplace_unique(v, (increment)); \
        } \
        \
        if (!(is_head && !w_in_set)) { \
            auto it_w = local_triangles.find(w); \
            if (it_w != local_triangles.end()) it_w->second += (increment); \
            else local_triangles.emplace_unique(w, (increment)); \
        } \
    } while(0)



    // Process main neighbors of u
    if (u_neighs) {
        for (const auto w: *u_neighs) {
            // Check if w is in v's main sample
            if (v_neighs && v_neighs->find(w) != v_neighs->end()) {
                UPDATE_TRIANGLES(w, inc_mm);
            }
                // Check if w is in v's aux sample
            else if (v_neighs_aux && v_neighs_aux->find(w) != v_neighs_aux->end()) {
                UPDATE_TRIANGLES(w, inc_ma);
            }
        }
    }

    // Process aux neighbors of u
    if (u_neighs_aux) {
        for (const auto w: *u_neighs_aux) {
            // Check if w is in v's main sample
            if (v_neighs && v_neighs->find(w) != v_neighs->end()) {
                UPDATE_TRIANGLES(w, inc_ma);
            }
                // Check if w is in v's aux sample
            else if (v_neighs_aux && v_neighs_aux->find(w) != v_neighs_aux->end()) {
                UPDATE_TRIANGLES(w, inc_aa);
            }
        }
    }

#undef UPDATE_TRIANGLES
    return cum_count;
}


/**
 * Sample edge {u, v} in head or tail or aux reservoir samplers
 * @param u
 * @param v
 */
void TriangleCounter::sample_edge(const int u, const int v) {
    const bool u_in_h = head_nodes_.find(u) != head_nodes_.end();
    const bool u_in_t = tail_nodes_.find(u) != tail_nodes_.end();
    const bool v_in_h = head_nodes_.find(v) != head_nodes_.end();
    const bool v_in_t = tail_nodes_.find(v) != tail_nodes_.end();

    // tail sampling
    if (u_in_t || v_in_t)
        reservoir_sample_edge(u, v, tail_sample_, sample_tail_cur_,
                              sample_tail_size_, tail_rand_indexer_, subgraph_tail_);
    else
        reservoir_sample_edge(u, v, aux_sample_tail_, aux_sample_cur_tail_,
                              aux_sample_size_tail_, aux_tail_rand_indexer_, aux_subgraph_tail_);


    // head sampling
    if (u_in_h || v_in_h)
        reservoir_sample_edge(u, v, head_sample_, sample_head_cur_,
                              sample_head_size_, head_rand_indexer_, subgraph_head_);
    else
        reservoir_sample_edge(u, v, aux_sample_head_, aux_sample_cur_head_,
                              aux_sample_size_head_, aux_head_rand_indexer_, aux_subgraph_head_);

}


/**
 * Process incoming edge {u, v}
 * @param u
 * @param v
 */
void TriangleCounter::process_edge(int u, int v) {

    if (u > v) std::swap(u, v);

    global_triangle_cnt_head_ += count_triangles(u, v, true, head_nodes_, subgraph_head_, aux_subgraph_head_,
                                                 sample_head_cur_, sample_head_size_, aux_sample_cur_head_,
                                                 aux_sample_size_head_,
                                                 local_triangle_head_cnt_);
    global_triangle_cnt_tail_ += count_triangles(u, v, false, tail_nodes_, subgraph_tail_, aux_subgraph_tail_,
                                                 sample_tail_cur_, sample_tail_size_, aux_sample_cur_tail_,
                                                 aux_sample_size_tail_,
                                                 local_triangle_tail_cnt_);
    sample_edge(u, v);


}