//
// Created by Cristian Boldrin on 14/11/25.
//

#include "Utils.h"

/**
 * Function that preprocesses a graph and saves it in the correct format, i.e., (u v t) for each row, separated
 * by a space delimiter. Also, handles self-loops and multiple edges and sorts the edges by increasing time of arrival
 * in the stream.
 * @param dataset_filepath path for the graph dataset file
 * @param delimiter for rows of graph dataset file
 * @param skip line to skip at the beginning of graph dataset file
 * @param output_path where to store the preprocess graph dataset
 */
void Utils::preprocess_data(const std::string &dataset_filepath, std::string &delimiter, int skip,
                            std::string &output_path) {

    std::cout << "Preprocessing Dataset...\n";
    std::ifstream file(dataset_filepath);
    std::string line, su, sv;

    // -- edge stream
    std::unordered_map<Edge, int, hash_edge> edge_stream;

    // - graph
    std::unordered_map<int, std::unordered_set<int>> graph_stream;

    int u, v, t;
    std::unordered_set<int> min_neighbors;

    if (file.is_open()) {

        long nline = 0;
        long num_nodes;
        long num_edges = 0;

        t = 0;
        while (std::getline(file, line)) {
            nline++;
            if (nline <= skip) continue;

            std::istringstream iss(line);
            std::getline(iss, su, delimiter[0]);
            std::getline(iss, sv, delimiter[0]);

            u = stoi(su);
            v = stoi(sv);

            // -- check self-loops
            if (u == v) continue;
            t++;
            int v1 = (u < v) ? u : v;
            int v2 = (u < v) ? v : u;
            std::pair uv = std::make_pair(v1, v2);
            // -- check for multiple edges
            if (graph_stream[u].find(v) != graph_stream[u].end() and
                graph_stream[v].find(u) != graph_stream[v].end()) {
                edge_stream[uv] = t;
                continue;
            }

            // -- add edge to graph stream
            graph_stream[u].emplace(v);
            graph_stream[v].emplace(u);
            num_edges++;
            edge_stream[uv] = t;

            if (nline % 3000000 == 0) {
                std::cout << "Processed " << nline << " edges...\n";
            }

        }

        // -- eof
        num_nodes = (int) graph_stream.size();
        printf("Preprocessed dataset with n = %ld, m = %ld\n", num_nodes, num_edges);
        std::cout << "Sorting edge map...\n";
        // -- create a vector that stores all the entries <K, V> of the map edge stream
        std::vector<std::pair<Edge, int>> ordered_edge_stream(edge_stream.begin(), edge_stream.end());
        // -- sort edge by increasing time
        std::sort(ordered_edge_stream.begin(), ordered_edge_stream.end(),
                  [](const std::pair<Edge, int> &a, const std::pair<Edge, int> &b) { return a.second < b.second; });

        // -- write results
        std::cout << "Done!\nWriting results...\n";
        std::ofstream out_file(output_path);

        int cnt = 0;
        for (auto elem: ordered_edge_stream) {
            // -- also, rescale the time (not meant for Tonic)
            out_file << elem.first.first << " " << elem.first.second << " " << ++cnt << "\n";
        }

        out_file.close();

    } else {
        std::cerr << "DataPreprocessing - Error! Graph filepath not opened.\n";
    }

}

/**
 * Function that runs the exact algorithm to count triangles and degrees in a single pass over the graph stream.
 * @param graph_path: path to the input graph stream (edge list format)
 * @param output_path: path to the output file (it will be a csv file with columns: node_id, degree, triangles)
 */
void Utils::run_exact_stream(const std::string &graph_path, const std::string &output_path) {

    GraphStream stream(graph_path);
    GraphStream::Edge current_edge{};
    emhash5::HashMap<int, emhash8::HashSet<int>> graph;
    emhash5::HashMap<int, int> triangle_counter;

    // -- main for loop: process edges
    while(stream.next(current_edge)) {

        size_t edge_count = stream.getEdgeCount();
        if (edge_count % 5000000 == 0)
            printf("Processed %zuM edges...\n", edge_count / 1000000);

        int u = (int) current_edge.u;
        int v = (int) current_edge.v;

        // -- lambda function to account node-triangle
        auto UPDATE_TRIANGLES = [&](int w) {
            // -- update counters local triangles
            auto w_it = triangle_counter.find(w);
            if (w_it != triangle_counter.end()) w_it->second += 1;
            else triangle_counter.emplace(w, 1);
        };

        // -- add to graph
        graph[u].emplace(v);
        graph[v].emplace(u);
        // -- count triangles
        auto du = (int) graph[u].size();
        auto dv = (int) graph[v].size();
        int n_min = (du < dv) ? u : v;
        int n_max = (du < dv) ? v : u;
        auto &min_neighbors = graph[n_min];
        for (const auto &neigh: min_neighbors) {
            if (graph[n_max].find(neigh) != graph[n_max].end()) {
                // -- triangle {n_min, neigh, n_max} discovered
                UPDATE_TRIANGLES(n_min);
                UPDATE_TRIANGLES(n_max);
                UPDATE_TRIANGLES(neigh);
            }
        }

    }

    // -- write estimates at the end
    std::ofstream out_file(output_path);
    // -- insert header
    out_file << "node_id,degree,triangles\n";
    for (const auto &node_it: graph) {
        auto deg = (int) node_it.second.size();
        auto tri_it = triangle_counter.find(node_it.first);
        int tri_count = (tri_it != triangle_counter.end()) ? tri_it->second : 0;
        out_file << node_it.first << "," << deg << "," << tri_count << "\n";
    }
    out_file.close();

}