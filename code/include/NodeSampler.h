//
// Created by Cristian Boldrin on 09/12/25.
//

#pragma once

#include <iostream>
#include <string>
#include <utility>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <vector>

#include <utils/hash_table5.hpp>
#include <utils/hash_set8.hpp>
#include <utils/random.h>

using node_degree_t = emhash5::HashMap<int, int>; // node -> degree

class NodeSampler {

private:

    double ph_, pt_, eps_;
    size_t random_seed_;

    // -- sample and counters
    node_degree_t head_sampled_nodes_, tail_sampled_nodes_; // node -> degree

    // -- deterministic seed for head hashing
    std::uint64_t head_hash_seed_ = 0;

    // -- expected degrees cache
    emhash5::HashMap<int, int> exp_deg_cache_;

    // -- rng
    utils::Random rand_gen_;

    static inline std::uint64_t splitmix64(std::uint64_t x) noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    /**
     * Deterministic hash function for head node sampling.
     * Returns a value in [0, 1) without storing anything per node.
     */
    inline double hasher(int node) const noexcept {
        const std::uint64_t x =
                (static_cast<std::uint64_t>(static_cast<std::uint32_t>(node)) << 32) ^
                static_cast<std::uint64_t>(static_cast<std::uint32_t>(node));
        const std::uint64_t z = splitmix64(x ^ head_hash_seed_);
        return static_cast<double>(z >> 11) * (1.0 / 9007199254740992.0); // 2^53
    }

    double expcount(double deg, double p) const;

    double expdegBS(double c, double p, double low, double high) const;

    int expdeg(int c, double p);

    int expdeg_const(int c, double p) const;

    emhash5::HashMap<int, int> compute_histogram(const emhash5::HashMap<int, int> &map) const;

public:

    /**
     * Sampling result for an edge
     */
    struct EdgeSamplingResult {
        bool u_sampled_head;
        bool v_sampled_head;
        bool u_sampled_tail;
        bool v_sampled_tail;
        bool u_newly_added_head;  // Was u added to head in this call?
        bool v_newly_added_head;
        bool u_newly_added_tail;
        bool v_newly_added_tail;
    };

    NodeSampler(double ph, double pt, double eps, size_t random_seed);

    ~NodeSampler();

    EdgeSamplingResult process_edge(const int u, const int v);

    void correct_tail_estimates();


    int compute_thresh() const;

    void print_stats() {
        std::cout << "[NodeSampler] Head sampled nodes |S_h| = " << head_sampled_nodes_.size() << "\n";
        std::cout << "[NodeSampler] Tail sampled nodes |S_t| = " << tail_sampled_nodes_.size() << "\n";
    }

    void get_head_sampled_nodes(node_degree_t &head_map) const {
        for (const auto &it : head_sampled_nodes_)
            head_map.insert_unique(it.first, it.second);
    }

    void get_tail_sampled_nodes(node_degree_t &tail_map) const {
        for (const auto &it : tail_sampled_nodes_)
            tail_map.insert_unique(it.first, it.second);
    }
};