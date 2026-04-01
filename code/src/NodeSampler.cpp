//
// Created by Cristian Boldrin on 09/12/25.
//

#include "NodeSampler.h"

NodeSampler::NodeSampler(double ph, double pt, double eps, size_t random_seed) :
        ph_(ph),
        pt_(pt),
        eps_(eps),
        random_seed_(random_seed),
        rand_gen_(random_seed) {

    // deterministic per-node hashing seed for head sampling
    head_hash_seed_ = splitmix64(static_cast<std::uint64_t>(random_seed_));

    printf("\n----- Initializing NodeSampler (headtail algorithm)... -----\n");
    printf("Head (Nodes) Sampling Probability p_h: %.4f\n", ph_);
    printf("Tail (Nodes) Sampling Probability p_t: %.4f\n", pt_);
}

NodeSampler::~NodeSampler() {
    head_sampled_nodes_.clear();
    tail_sampled_nodes_.clear();
    exp_deg_cache_.clear();
}

/**
 * The count computed as d(v) - E[loss(v) | v is sampled]
 *
 * c = d*p + p - 1 + (1-p)^(d+1) / (p * (1 - (1-p)^d))
 * @param deg
 * @param p
 * @return
 */
double NodeSampler::expcount(double deg, double p) const {
    auto num = deg * p + p - 1 + std::pow((1 - p), (deg + 1));
    auto den = p * (1 - std::pow((1 - p), deg));
    return num / den;
}

/**
 * Use binary search to find degree d such that expcount(d, p) = c
 * @param c
 * @param p
 * @param low
 * @param high
 * @return
 */
double NodeSampler::expdegBS(double c, double p, double low, double high) const {
    double guess = (low + high) / 2.0;
    double deg_val = expcount(guess, p);
    auto rounded_deg_val = std::round(deg_val);
    auto rounded_c = std::round(c);

    if (rounded_deg_val == rounded_c) {
        return guess;
    }
    else if (rounded_deg_val < rounded_c) return expdegBS(c, p, guess, high);
    else if (rounded_deg_val > rounded_c) return expdegBS(c, p, low, guess);
    else {
        std::cerr << "Error in expdegBS()\n";
        return -1.0;
    }
}

/**
 * Invert expcount to get expected degree from expected count
 * @param c
 * @param p
 * @return
 */
int NodeSampler::expdeg(int c, double p) {
    auto cache_it = exp_deg_cache_.find(c);
    if (cache_it != exp_deg_cache_.end()) return cache_it->second;

    auto deg = static_cast<int>(std::round(expdegBS(static_cast<double>(c), p, c, 2 * c)));
    exp_deg_cache_.emplace_unique(c, deg);
    return deg;
}

/**
 * Same as expdeg but without caching (const version)
 * @param c
 * @param p
 * @return
 */
int NodeSampler::expdeg_const(int c, double p) const {
    auto deg = static_cast<int>(std::round(expdegBS(static_cast<double>(c), p, c, 2 * c)));
    return deg;
}

/**
 * Correct tail node degree estimates using expdeg function
 */
void NodeSampler::correct_tail_estimates() {
    for (auto &it : tail_sampled_nodes_) {
        int sampled_count = it.second;
        int corrected_deg = expdeg(sampled_count, pt_);
        it.second = corrected_deg;
    }
}


emhash5::HashMap<int, int> NodeSampler::compute_histogram(const emhash5::HashMap<int, int> &map) const {
    emhash5::HashMap<int, int> histogram;
    for (const auto &it : map) {
        int deg = it.second;
        auto hist_it = histogram.find(deg);
        if (hist_it != histogram.end()) hist_it->second += 1;
        else histogram.emplace_unique(deg, 1);
    }
    return histogram;
}

int NodeSampler::compute_thresh() const {
    auto head_hist = compute_histogram(head_sampled_nodes_);

    std::vector<std::pair<int, int>> head_hist_vec;
    for (const auto &it : head_hist)
        head_hist_vec.emplace_back(it.first, it.second);

    std::sort(head_hist_vec.begin(), head_hist_vec.end());

    auto thresh = (3.0 * std::log(1.0 / eps_)) / (ph_ * std::pow(eps_, 2));
    int deg_thresh = 0;

    for (const auto& it : head_hist_vec) {
        int deg = it.first;
        int count = it.second;
        double nd = count / ph_;
        if (nd >= thresh) deg_thresh = deg;
    }

    return deg_thresh;
}

NodeSampler::EdgeSamplingResult NodeSampler::process_edge(const int u, const int v) {

    EdgeSamplingResult result = {false, false, false, false, false, false, false, false};

    // -- head
    auto u_head_it = head_sampled_nodes_.find(u);
    auto v_head_it = head_sampled_nodes_.find(v);

    if (u_head_it != head_sampled_nodes_.end()) {
        u_head_it->second += 1;
        result.u_sampled_head = true;
    }
    else {
        double alpha_1 = hasher(u);
        if (alpha_1 < ph_) {
            head_sampled_nodes_.emplace_unique(u, 1);
            result.u_newly_added_head = true;
        }
    }

    if (v_head_it != head_sampled_nodes_.end()) {
        v_head_it->second += 1;
        result.v_sampled_head = true;
    }
    else {
        double beta_1 = hasher(v);
        if (beta_1 < ph_) {
            head_sampled_nodes_.emplace_unique(v, 1);
            result.v_newly_added_head = true;
        }
    }

    // -- tail
    auto u_tail_it = tail_sampled_nodes_.find(u);
    auto v_tail_it = tail_sampled_nodes_.find(v);

    if (u_tail_it != tail_sampled_nodes_.end()) {
        u_tail_it->second += 1;
        result.u_sampled_tail = true;
    }
    else {
        double alpha_2 = rand_gen_.getDouble();
        if (alpha_2 < pt_) {
            tail_sampled_nodes_.emplace_unique(u, 1);
            result.u_newly_added_tail = true;
        }
    }

    if (v_tail_it != tail_sampled_nodes_.end()) {
        v_tail_it->second += 1;
        result.v_sampled_tail = true;
    }
    else {
        double beta_2 = rand_gen_.getDouble();
        if (beta_2 < pt_) {
            tail_sampled_nodes_.emplace_unique(v, 1);
            result.v_newly_added_tail = true;
        }
    }

    return result;
}