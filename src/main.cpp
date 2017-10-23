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

std::map<std::string, std::map<std::string, int> > load_mappings(std::string mapping_filename) {
    std::ifstream infile(mapping_filename);
    std::string line;
    std::string fields[3];
    std::string barcode;
    int counter = 0;
    std::cout << "Loading mappings file " << mapping_filename << std::endl;
    std::map<std::string, std::map<std::string, int> > mappings;
    while (std::getline(infile, line)){
        // read name, contig, number of kmers
        std::istringstream(line) >> fields[0] >> fields[1] >> fields[2] ;
        barcode = fields[0].substr(fields[0].find("_") + 1);
        mappings[fields[2]][barcode] += std::stoi(fields[1]);
        counter += 1;

    }
    std::cout << "Loaded " << counter << " mappings  " << std::endl;
    return mappings;
}

int main(int argc, char **argv) {

    const unsigned int GB(1024 * 1024 * 1024);
    unsigned int mem_limit(10);

    std::string graph_file_list;
    std::string mappings_filename;
    std::string output_file_pref;

    graph_file_list = argv[1];
    mappings_filename = argv[2];
    output_file_pref = argv[3];
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
    // actually should do it with choice of single gfa or mapping file list
    // loading entire mappings file each time takes ages- better strategy is to lload whole thing
    // take list of gfas, and loop from in here
    std::map<std::string, std::map<std::string, int> > mappings = load_mappings(mappings_filename);
    std::string graph_filename;
    std::string start_edge;
    std::string fields[2];
    std::string line;
    std::ifstream infile(graph_file_list);
    // this should be a oarallel for
    # pragma omp parallel for
    while (std::getline(infile, line)) {
        std::istringstream(line) >> fields[0] >> fields[1];
        graph_filename = fields[0];
        start_edge = fields[1];
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "Phasing GFA: " << graph_filename << std::endl;
        std::string filename = graph_filename.substr(2, graph_filename.find_last_of(".gfa"));
        std::string output_file = output_file_pref + filename;
        std::cout << "Output file: " << output_file <<std::endl;
        Graph graph = Graph();
        graph.load_gfa(graph_filename);
        std::cout << "Traversing from start edge " << start_edge << " in + direction" << std::endl;
        std::vector<std::string> traversed_edge_list;
        graph.traverse_graph(start_edge, "+", traversed_edge_list);
        std::cout << "Found " << graph.bubbles.size() << " bubbles from + direction" << std::endl;
        std::cout << "Traversing from start edge " << start_edge << " in - direction" << std::endl;
        traversed_edge_list.clear();
        graph.traverse_graph(start_edge, "-", traversed_edge_list);

        std::cout << "Found " << graph.bubbles.size() << " bubbles  in total" << std::endl;
        if (graph.bubbles.size() > 1) {
            std::vector<std::vector<std::string> > possible_haplotypes = graph.calculate_possible_haplotypes();
            std::cout << "found " << possible_haplotypes.size() << "candidate haplotypes of length "
                      << possible_haplotypes[0].size() << std::endl;
            std::cout << "loading " << mappings_filename << " " << mappings.size() << std::endl;
            HaplotypeScorer haplotype_scorer = HaplotypeScorer(mappings_filename, possible_haplotypes, graph);
            haplotype_scorer.load_mappings_from_dict(mappings);
            haplotype_scorer.decide_barcode_haplotype_support();
            int success = haplotype_scorer.score_haplotypes("formatted_" + output_file);

            // if we've picked a winner
            if (success == 0) {
                std::cout << "Writing output" <<std::endl;
                haplotype_scorer.write_output_success(output_file);
                graph.write_output_subgraph(std::get<0>(haplotype_scorer.winners), "sequences1" + output_file + ".gfa",
                                            "haplotype1");
                graph.write_output_subgraph(std::get<1>(haplotype_scorer.winners), "sequences2" + output_file + ".gfa",
                                            "haplotype2");

            } else if (success == 1) { // if we're less confident about winner
                std::cout << "Writing output" <<std::endl;

                haplotype_scorer.write_output_partial_success(output_file);
                graph.write_output_subgraph(std::get<0>(haplotype_scorer.winners),
                                            "partial_sequences1" + output_file + ".gfa",
                                            "haplotype1");
                graph.write_output_subgraph(std::get<1>(haplotype_scorer.winners),
                                            "partial_sequences2" + output_file + ".gfa",
                                            "haplotype2");

            }
        }
        std::cout << "----------------------------------------" << std::endl;
        std::cout <<std::endl;
        std::cout <<std::endl;
        std::cout <<std::endl;
    }
    return 0;
}
