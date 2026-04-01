#include <iostream>
#include <string>
#include <Utils.h>
#include <GraphStream.h>
#include <TriangleCounter.h>
#include <NodeSampler.h>

/**
 * Helper function to write node estimates (degrees or triangles) to file in csv format (node_id, estimate)
 * @tparam ValueType: type of the estimate (e.g., int for degrees, double for triangle counts)
 * @param output_path: path to the output file
 * @param overwrite: whether to overwrite the file (true) or append to it (false)
 * @param node_estimates: hash map containing node estimates (node_id -> estimate)
 */
template<typename ValueType>
void write_results_to_file(const std::string &output_path, bool overwrite,
                           const emhash5::HashMap<int, ValueType> &node_estimates) {
    std::ofstream out_file;
    if (overwrite) out_file.open(output_path); // overwrite mode
    else out_file.open(output_path, std::ios::app); // append mode
    for (const auto &it: node_estimates)
        out_file << it.first << "," << it.second << "\n";
    out_file.close();
}

/**
 * Get the base name of the executable
 * @param s the string to split
 * @return the name of the executable
 */
char *base_name(char *s) {
    char *start;
    if ((start = strrchr(s, '/')) == nullptr) start = s;
    else ++start;
    return start;
}

int main(int argc, char **argv) {

    char *project = base_name(argv[0]);

    // -- data preprocessing
    if (strcmp(project, "DataPreprocessing") == 0) {
        if (argc != 5) {
            printf("Invalid number of arguments.\n");
            printf("Usage: DataPreprocessing <dataset_path> <delimiter> <skip> <output_path>\n");
            std::cout << "\t- dataset_path: path to the input dataset (edge list format)\n";
            std::cout << "\t- delimiter: delimiter used in the input dataset (e.g., space, comma, tab), e.g., for tab, use $'\\t'\n";
            std::cout << "\t- skip: number of lines to skip at the beginning of the file (e.g., for header)\n";
            std::cout << "\t- output_path: path to the output file\n";
            return 1;
        } else {
            std::string dataset_path(argv[1]);
            std::string delimiter (argv[2]);
            int skip = atoi(argv[3]);
            std::string output_path(argv[4]);
            utils::StopWatch sw;
            sw.start();
            Utils::preprocess_data(dataset_path, delimiter, skip, output_path);
            auto time = sw.elapsedSeconds();
            std::cout << "Dataset preprocessed in time: " << time << " s\n";
            return 0;
        }
    }

    // -- Run Exact (Count Triangles and Degrees exactly in a single pass)
    if (strcmp(project, "RunExact") == 0) {
        if (argc != 3) {
            printf("Invalid number of arguments.\n");
            printf("Usage: RunExact <graph_path> <output_path>\n");
            std::cout << "\t- graph_path: path to the input graph stream (edge list format)\n";
            std::cout << "\t- output_path: path to the output file (it will be a csv file)\n";
            return 1;
        } else {
            std::string graph_path(argv[1]);
            std::string output_path(argv[2]);
            std::cout << "Running exact algorithm on stream...\n";
            utils::StopWatch sw;
            sw.start();
            Utils::run_exact_stream(graph_path, output_path);
            auto time = sw.elapsedSeconds();
            std::cout << "Exact algorithm completed in time: " << time << " s\n";
            return 0;
        }
    }

    if (strcmp(project, "BOLIDE") == 0) {

        if (argc != 11) {
            printf("Invalid number of arguments.\n");
            printf("Usage: BOLIDE <graph_path> <p_head> <p_tail> <eps> "
                   "<output_path> <random_seed_sample> <random_seed_triangle_counter> "
                   "<aux_sample_size> <head_budget> <tail_budget>\n");
            std::cout << "\t- graph_path: path to the input graph stream (edge list format)\n";
            std::cout << "\t- p_head: head sampling probability (for nodes) via hashing\n";
            std::cout << "\t- p_tail: tail sampling probability (for nodes) via sample-and-hold\n";
            std::cout << "\t- eps: error parameter for degree threshold. From headtail, \\tau = 3.0*log(1/eps) / (eps**2); "
                         "In the paper, we used \\tau = 10, thus eps ~ 0.5\n";
            std::cout << "\t- output_path: path to the output file (it will be a csv file)\n";
            std::cout << "\t- random_seed_sample: random seed for sampling nodes (head and tail)\n";
            std::cout << "\t- random_seed_triangle_counter: random seed for triangle counter (sampling edges)\n";
            std::cout << "\t- aux_sample_size: size of the auxiliary edge sample (to be split in half for head and tail)\n";
            std::cout << "\t- head_budget: budget for the head edge sample (number of edges)\n";
            std::cout << "\t- tail_budget: budget for the tail edge sample (number of edges)\n";

            return 1;
        }

        // -- read input parameters
        std::string graph_path(argv[1]);
        double ph = std::stod(argv[2]);
        double pt = std::stod(argv[3]);
        double eps = std::stod(argv[4]);
        std::string output_path(argv[5]);
        size_t random_seed_sampler = std::stol(argv[6]);
        size_t random_seed_tc = std::stol(argv[7]);
        size_t aux_sample_size = std::stol(argv[8]);
        size_t head_budget = std::stol(argv[9]);
        size_t tail_budget = std::stol(argv[10]);

        NodeSampler sampler(ph, pt, eps, random_seed_sampler);
        // -- aux sample size is split in half for head and tail
        // -- (design choice, we could also have separate parameters for head and tail aux sample sizes)
        TriangleCounter triangle_counter(aux_sample_size / 2,
                                         aux_sample_size - aux_sample_size / 2,
                                         head_budget, tail_budget,random_seed_tc);

        std::string line;
        utils::StopWatch sw;
        sw.start();

        // -- initialize buffer for input stream
        GraphStream stream(graph_path);
        GraphStream::Edge current_edge{};

        // -- main for loop: process each incoming edge in the stream
        while (stream.next(current_edge)) {

            size_t edge_count = stream.getEdgeCount();
            // -- log progress every 1M edges
            if (edge_count % 1000000 == 0)
                printf("Processed %zuM edges...\n", edge_count / 1000000);


            int u = (int) current_edge.u;
            int v = (int) current_edge.v;

            auto sample_result = sampler.process_edge(u, v);
            if (sample_result.u_newly_added_head) triangle_counter.add_node_head_sample(u);
            if (sample_result.v_newly_added_head) triangle_counter.add_node_head_sample(v);
            if (sample_result.u_newly_added_tail) triangle_counter.add_node_tail_sample(u);
            if (sample_result.v_newly_added_tail) triangle_counter.add_node_tail_sample(v);
            triangle_counter.process_edge(u, v);

        }

        // -- correct tail estimates (correct degree estimates using sample-and-hold correction, i.e., expected loss of TG distribution)
        sampler.correct_tail_estimates();
        auto deg_thresh = sampler.compute_thresh();

        auto time = sw.elapsedSeconds();
        std::cout << "Done!\n-----> Time elapsed (seconds): " << time << "\n";

        // -- log useful statistics about the samples
        sampler.print_stats();

        // -- fetch estimates, and log stats
        emhash5::HashMap<int, double> triangles_head, triangles_tail;
        triangle_counter.get_estimates_head(triangles_head, true);
        triangle_counter.get_estimates_tail(triangles_tail, true);
        auto sample_size = triangle_counter.get_subgraph_size();
        std::cout << "-----> Final sample size: " << sample_size << " edges.\n";
        auto unique_edges = triangle_counter.get_unique_edges();
        std::cout << "-----> Unique edges in sample: " << unique_edges << " edges.\n";

        // -- write to file
        std::string stem = output_path.substr(0, output_path.find_last_of('.'));
        // -- results sampler
        emhash5::HashMap<int, int> head_sampled_nodes, tail_sampled_nodes;
        sampler.get_head_sampled_nodes(head_sampled_nodes);
        sampler.get_tail_sampled_nodes(tail_sampled_nodes);
        // -- write node-degrees to file
        // -- head
        std::string head_deg_output_path = stem + "_head_node_degrees.txt";
        write_results_to_file(head_deg_output_path, true, head_sampled_nodes);
        // -- tail
        std::string tail_deg_output_path = stem + "_tail_node_degrees.txt";
        write_results_to_file(tail_deg_output_path, true, tail_sampled_nodes);


        // -- results triangle counter
        // -- head
        std::string head_output_path = stem + "_head_node_triangles.txt";
        write_results_to_file(head_output_path, true, triangles_head);
        // -- tail
        std::string tail_output_path = stem + "_tail_node_triangles.txt";
        write_results_to_file(tail_output_path, true, triangles_tail);


        // -- save some additional infos
        std::string info_path = stem + "_info.csv";
        std::ofstream info_out_file(info_path, std::ios::app);
        // -- sample_size,unique_sample_size,aux_sample_size,head_budget,tail_budget,degree_threshold,time
        info_out_file << sample_size << "," << unique_edges << "," << aux_sample_size << "," << head_budget << ","
                      << tail_budget << "," << deg_thresh << "," << time << "\n";
        info_out_file.close();

        return 0;


    }

}