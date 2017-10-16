#include <iostream>
#include <fstream>

#include <string>
#include <unordered_map>

#include <sys/stat.h>

//#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

#include "../deps/cxxopts/include/cxxopts.hpp"
#include "graph.h"
#include "haplotype_scorer.h"

Graph load_subgraph(std::string subgraph_filename){
    std::ifstream subgraph(subgraph_filename);

}

bool check_or_create_directory(std::string &output_prefix) {
    if (output_prefix.back() != '/') {
        output_prefix.push_back('/');
    }
    struct stat sb{};
    bool validate_dir(false);
    if (stat(output_prefix.c_str(), &sb) != 0) {
        if (errno == ENOENT) {
            mode_t mask = umask(0);
            umask(mask);
            mkdir(output_prefix.c_str(), mode_t(0777 - mask));
            validate_dir = true;
        }
        if (stat(output_prefix.c_str(), &sb) != 0) {
            perror(output_prefix.c_str());
            validate_dir = false;
        }
    } else if (!S_ISDIR(sb.st_mode)) {
        std::cout << output_prefix << " is not a directory " << std::endl;
    } else {
        validate_dir = true;
    }
    return validate_dir;
}

int main(int argc, char **argv) {

    const unsigned int GB(1024 * 1024 * 1024);
    unsigned int mem_limit(10);

    std::string graph_filename;
    std::string mappings_filename;
    std::string start_edge;
    std::string output_file;

    graph_filename = argv[1];
    start_edge = argv[2];
    mappings_filename = argv[3];
    /*try {
        std::time_t t = std::time(nullptr);
        std::tm tm = *std::localtime(&t);
        std::string outdefault(
                std::to_string(tm.tm_year + 1900) + '-' + std::to_string(tm.tm_mon) + '-' + std::to_string(tm.tm_mday) +
                '_' + std::to_string(tm.tm_hour) + std::to_string(tm.tm_min));

        cxxopts::Options options("phaser", "Phase part of a graph based on kmer mapping");

        options.add_options()("help", "Print help")("g, graph", "Graph to phase",
                                                    cxxopts::value<std::string>(graph_filename), "GFA - Graph file")(
                "m,mappings", "Mappings used to phase graph",
                cxxopts::value<std::string>(mappings_filename), "Mappings")("s,start",
                                                                            "Edge to start on",
                                                                            cxxopts::value<std::string>(
                                                                                    start_edge),
                                                                            "Start");

        options.add_options("Output file")("o,output", "Output proposed phasing",
                                           cxxopts::value<std::string>(output_file), "Output file");

        options.parse(argc, argv);

        if (0 != options.count("help")) {
            std::cout << options.help({"", "Performance"}) << std::endl;
            exit(0);
        }
        if (options.count("o") != 1 /*or options.count("i")<2) {
            std::cout << "Error: please specify input files and output prefix" << std::endl
                      << " Use option --help to check command line arguments." << std::endl;
            exit(1);
        }


        if (graph_filename.empty()) {
            std::cout << "Error: The GFA file parameter wasn't specified, " << std::endl
                      << "Use option --help to check command line arguments." << std::endl;
            exit(1);
        }
        if (mappings_filename.empty()) {
            std::cout << "Error: The mapping file parameter wasn't specified, " << std::endl
                      << "Use option --help to check command line arguments." << std::endl;
            exit(1);
        }
    } catch (const cxxopts::OptionException &e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }*/

    Graph graph=Graph();
    graph.load_gfa(graph_filename);
    std::cout << "Traversing from start edge " << start_edge << " in + direction" << std::endl;
    std::vector<std::string > traversed_edge_list;
    graph.traverse_graph(start_edge, "+", traversed_edge_list);
    std::cout << "Found " << graph.bubbles.size() << " bubbles from + direction" << std::endl;
    std::cout << "Traversing from start edge " << start_edge << " in - direction" << std::endl;
    traversed_edge_list.clear();
    graph.traverse_graph(start_edge, "-", traversed_edge_list);

    std::cout << "Found " << graph.bubbles.size() << " bubbles  in total" << std::endl;
    for (auto l:graph.bubbles){
        std::cout << std::get<0>(l) << " " << std::get<1>(l) << std::endl;
    }
    std::vector<std::vector <std::string> > possible_haplotypes = graph.calculate_possible_haplotypes();
    std::cout<< "found " << possible_haplotypes.size() << "candidate haplotypes of length " << possible_haplotypes[0].size() << std::endl;
    for (auto h: possible_haplotypes){
        for (auto i:h){
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
    std::cout<<  "loading " << mappings_filename << std::endl;
    HaplotypeScorer haplotype_scorer = HaplotypeScorer(mappings_filename, possible_haplotypes, graph);
    haplotype_scorer.load_mappings();
    haplotype_scorer.decide_barcode_haplotype_support();
    std::pair<std::vector<std::string>,std::vector<std::string> >  winners  = haplotype_scorer.score_haplotypes();
    return 0;
}