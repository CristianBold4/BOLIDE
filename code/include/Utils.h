//
// Created by Cristian Boldrin on 14/11/25.
//

#pragma once

#include <iostream>
#include <string>
#include <utils/stop_watch.h>
#include <utils/random.h>

#include "utils/hash_table5.hpp"
#include "utils/hash_set8.hpp"
#include "GraphStream.h"

class Utils {

    using Edge = std::pair<int, int>;

    constexpr static unsigned long long MAX_ID_NODE = 100000000;

    inline static unsigned long long edge_to_id(const Edge e) {
        int u = e.first;
        int v = e.second;
        int nu = (u < v ? u : v);
        int nv = (u < v ? v : u);
        return (MAX_ID_NODE) * static_cast<unsigned long long>(nu) +
               static_cast<unsigned long long>(nv);
    }

    struct hash_edge {
        size_t operator()(const Edge &e) const {
            return edge_to_id(e);
        }
    };

public:

    static void preprocess_data(const std::string &dataset_filepath, std::string &delimiter, int skip,
                                std::string &output_path);

    static void run_exact_stream(const std::string& graph_path, const std::string& output_path);

};


