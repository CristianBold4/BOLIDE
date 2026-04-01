//
// Created by Cristian Boldrin on 17/11/25.
//

#pragma once

#include <iostream>
#include <string>
#include <utility>
#include <utils/random.h>
#include <utils/hash_table5.hpp>
#include <utils/hash_set8.hpp>

using edge_t = std::pair<int, int>;
using graph_t = emhash5::HashMap<int, emhash8::HashSet<int>>;

using triangle_estimate_t = emhash5::HashMap<int, double>;
using node_triangle_t = triangle_estimate_t;
using my_hash_set_t = emhash8::HashSet<int>;

class TriangleCounter {

private:

    // -- sg: weight on edges indicates whether the edge is in a main sample, or in the aux sample
    graph_t subgraph_head_;
    graph_t subgraph_tail_;
    graph_t aux_subgraph_head_;
    graph_t aux_subgraph_tail_;

    // -- aux sample
    std::vector<edge_t> aux_sample_head_;
    std::vector<edge_t> aux_sample_tail_;
    size_t aux_sample_cur_head_ = 0;
    size_t aux_sample_cur_tail_ = 0;


    my_hash_set_t head_nodes_;
    my_hash_set_t tail_nodes_;

    utils::Random rand_gen_;
    utils::RandomIndexer aux_head_rand_indexer_;
    utils::RandomIndexer aux_tail_rand_indexer_;
    utils::RandomIndexer head_rand_indexer_;
    utils::RandomIndexer tail_rand_indexer_;

    double global_triangle_cnt_head_ = 0.0;
    double global_triangle_cnt_tail_ = 0.0;
    triangle_estimate_t local_triangle_head_cnt_;
    triangle_estimate_t local_triangle_tail_cnt_;

    size_t sample_head_cur_ = 0;
    size_t sample_tail_cur_ = 0;

    // -- vector for reservoir sampling (fast eviction)
    std::vector<edge_t> head_sample_;
    std::vector<edge_t> tail_sample_;

    struct EdgeHash {
        std::size_t operator()(edge_t e) const noexcept {
            if (e.first > e.second) std::swap(e.first, e.second);
            // Simple 64-bit mix
            std::uint64_t x = (static_cast<std::uint64_t>(e.first) << 32)
                              ^ static_cast<std::uint64_t>(e.second);
            x ^= x >> 33;
            x *= 0xff51afd7ed558ccdULL;
            x ^= x >> 33;
            x *= 0xc4ceb9fe1a85ec53ULL;
            x ^= x >> 33;
            return static_cast<std::size_t>(x);
        }
    };

    double count_triangles(int u, int v, bool is_head, const my_hash_set_t &set,
                           const graph_t &subgraph, const graph_t &aux_subgraph,
                           size_t cur, size_t size, size_t aux_cur, size_t aux_size,
                           node_triangle_t &local_triangles);

    void sample_edge(const int u, const int v);

    void add_edge(int u, int v, graph_t &graph);

    bool delete_edge(int u, int v, graph_t &graph);

    inline void reservoir_sample_edge(
            int u, int v,
            std::vector<std::pair<int, int>> &sample,
            size_t &cur,
            const size_t size,
            utils::RandomIndexer &rand_indexer,
            graph_t &subgraph
    ) {
        if (cur < size) {
            // Fill the reservoir
            sample[cur++] = {u, v};
            add_edge(u, v, subgraph);
        } else {
            cur++;
            // Reservoir sampling
            double p = (double) size / (double) cur;
            if (rand_gen_.getDouble() < p) {
                // Replace a random element
                int replace_idx = (int) rand_indexer.next();
                auto evicted_edge = sample[replace_idx];
                sample[replace_idx] = {u, v};
                // Remove evicted edge and add new one
                delete_edge(evicted_edge.first, evicted_edge.second, subgraph);
                add_edge(u, v, subgraph);
            }
        }
    }


public:

    size_t aux_sample_size_head_;
    size_t aux_sample_size_tail_;
    size_t sample_head_size_;
    size_t sample_tail_size_;

    TriangleCounter(size_t aux_head_sample_size, size_t aux_tail_sample_size,
                    size_t head_budget, size_t tail_budget, size_t random_seed);

    ~TriangleCounter();

    void process_edge(int u, int v);

    inline size_t get_subgraph_size() const {
        size_t sample_size = 0;
        for (const auto &it: subgraph_head_) sample_size += it.second.size();
        for (const auto &it: subgraph_tail_) sample_size += it.second.size();
        for (const auto &it: aux_subgraph_head_) sample_size += it.second.size();
        for (const auto &it: aux_subgraph_tail_) sample_size += it.second.size();
        return sample_size / 2; // -- each edge counted twice
    }

    inline size_t get_unique_edges() {
        std::unordered_set<edge_t, EdgeHash> unique_edges;
        auto process_graph = [&](const graph_t &graph) {
            for (const auto &it: graph) {
                int u = it.first;
                for (const auto &v: it.second) {
                    if (v <= u) continue; // -- ensure u < v (asymmetric checks)
                    unique_edges.emplace(u, v);
                }
            }
        };
        process_graph(subgraph_head_);
        process_graph(subgraph_tail_);
        process_graph(aux_subgraph_head_);
        process_graph(aux_subgraph_tail_);
        return unique_edges.size();

    }


    inline void get_estimates_head(node_triangle_t& estimates, bool verbose = false) const {
        double check_sum = 0.0;
        for (const auto &it: local_triangle_head_cnt_) {
            int node_id = it.first;
            double triangle_estimate = it.second;
            if (head_nodes_.find(node_id) == head_nodes_.end()) continue; // -- track only head nodes
            check_sum += triangle_estimate;
            estimates.emplace_unique(node_id, triangle_estimate);
        }

        if (!verbose) return;
        // -- capacity of head sample
        double capacity = std::min(1.0, (double) sample_head_cur_ / (double) sample_head_size_);
        // printf("> [HEAD] Check sum of local triangle estimates = %.3f\n", check_sum);
        printf("> [HEAD] Local Triangles estimates (||\\hat{T}_h||) = %u\n", local_triangle_head_cnt_.size());
        printf("> [HEAD] Size of estimates on nodes with tracked degree (S_h \\cup |\\hat{T}_h|) = %u\n", estimates.size());
        // printf("> [HEAD] Head sample capacity = %.3f%%\n", capacity * 100.0);
    }

    inline void get_estimates_tail(node_triangle_t &estimates, bool verbose = false) const {
        double check_sum = 0.0;
        for (const auto &it: local_triangle_tail_cnt_) {
            int node_id = it.first;
            double triangle_estimate = it.second;
            if (tail_nodes_.find(node_id) == tail_nodes_.end()) continue; // -- track only tail nodes
            check_sum += triangle_estimate;
            estimates.emplace_unique(node_id, triangle_estimate);
        }

        if (!verbose) return;
        double capacity = std::min(1.0, (double) sample_tail_cur_ / (double) sample_tail_size_);
        // printf("> [TAIL] Check sum of local triangle estimates = %.3f\n", check_sum);
        printf("> [TAIL] Local Triangles estimates (|\\hat{T}_t|) = %u\n", local_triangle_tail_cnt_.size());
        printf("> [TAIL] Size of estimates on nodes with tracked degree (S_t \\cup |\\hat{T}_t|) = %u\n", estimates.size());
        // printf("> [TAIL] Tail sample capacity = %.3f%%\n", capacity * 100.0);
    }


    void add_node_head_sample(const int u) { head_nodes_.insert(u);}

    void add_node_tail_sample(const int u) { tail_nodes_.insert(u);}


};

