#include <iostream>
#include <fstream>

#include <string>
#include <unordered_map>

#include <sys/stat.h>

//#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

#include "graph.h"
#include "haplotype_scorer.h"

bool check_or_create_directory(std::string &output_prefix) {
    if (output_prefix.empty()) {
        std::cout << "Error: output prefix must not be empty" << std::endl;
        return false;
    }
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
    // contig name -> barcode -> [number of kmers]
    std::map<std::string, std::map<std::string, int> > mappings;
    while (std::getline(infile, line)){
        // read name, number of kmerscontig name
        if (!(std::istringstream(line) >> fields[0] >> fields[1] >> fields[2])) {
            std::cout << "Warning: skipping malformed mappings line: " << line << std::endl;
            continue;
        }
        // extract barcode from read name, we know that reads with same barcode are from same chromosome
        barcode = fields[0].substr(fields[0].find("_") + 1);
        try {
            mappings[fields[2]][barcode] += std::stoi(fields[1]);
        } catch (const std::exception &e) {
            std::cout << "Warning: skipping mappings line with non-numeric kmer count: " << line << std::endl;
            continue;
        }
        counter += 1;

    }
    std::cout << "Loaded " << counter << " mappings  " << std::endl;
    return mappings;
}

int main(int argc, char **argv) {

    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " <graph_file_list> <mappings_file> <output_prefix>" << std::endl;
        return 1;
    }

    std::string graph_file_list(argv[1]);
    std::string mappings_filename(argv[2]);
    std::string output_file_pref(argv[3]);

    {
        std::ifstream check(mappings_filename);
        if (!check.is_open()) {
            std::cout << "Error: could not open mappings file " << mappings_filename << std::endl;
            return 1;
        }
    }
    {
        std::ifstream check(graph_file_list);
        if (!check.is_open()) {
            std::cout << "Error: could not open graph file list " << graph_file_list << std::endl;
            return 1;
        }
    }
    if (!check_or_create_directory(output_file_pref)) {
        std::cout << "Error: could not create or access output directory " << output_file_pref << std::endl;
        return 1;
    }

    // actually should do it with choice of single gfa or mapping file list
    // loading entire mappings file each time takes ages- better strategy is to lload whole thing
    // take list of gfas, and loop from in here
    std::map<std::string, std::map<std::string, int> > mappings = load_mappings(mappings_filename);
    std::string graph_filename;
    std::string start_edge;
    std::string fields[2];
    std::string line;
    std::ifstream infile(graph_file_list);
    int graphs = 0;
    int unphaseable = 0;
    int phaseable = 0;
    int phased_success = 0;
    int phased_partial_success = 0;
    int exceptions = 0;
    int no_mappings = 0;
    while (std::getline(infile, line)) {
        graphs += 1;
        if (!(std::istringstream(line) >> fields[0] >> fields[1])) {
            std::cout << "Warning: skipping malformed graph list line: " << line << std::endl;
            continue;
        }
        graph_filename = fields[0];
        start_edge = fields[1];
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "Phasing GFA: " << graph_filename << std::endl;
        int start = 0;
        if (graph_filename.find("/") != std::string::npos) {
            start = graph_filename.find("/");
        }
        std::string filename = graph_filename.substr(start + 1, graph_filename.find_last_of(".") - 1);
        std::string output_file = output_file_pref + filename;
        std::cout << "Output file: " << output_file << std::endl;
        Graph graph = Graph();
        graph.load_gfa(graph_filename);
        std::cout << "Traversing from start edge " << start_edge << " in + direction" << std::endl;
        std::vector<std::string> traversed_edge_list;
        // traverse forwards 
        graph.traverse_graph(start_edge, Strand::Plus, traversed_edge_list);
        std::cout << "Found " << graph.bubbles.size() << " bubbles from + direction" << std::endl;
        std::cout << "Traversing from start edge " << start_edge << " in - direction" << std::endl;
        traversed_edge_list.clear();
        // traverse backwards
        graph.traverse_graph(start_edge, Strand::Minus, traversed_edge_list);

        std::cout << "Found " << graph.bubbles.size() << " bubbles  in total" << std::endl;
        if (graph.bubbles.size() > 1) { // A length-1 haplotype can never produce more than 1 matching edge, so that condition can never be satisfied
            std::vector<std::vector<std::string> > possible_haplotypes = graph.calculate_possible_haplotypes();
            std::cout << "found " << possible_haplotypes.size() << "candidate haplotypes of length "
                      << possible_haplotypes[0].size() << std::endl;
            std::cout << "loading " << mappings_filename << " " << mappings.size() << std::endl;
            HaplotypeScorer haplotype_scorer = HaplotypeScorer(mappings_filename, possible_haplotypes, graph);
            haplotype_scorer.load_mappings_from_dict(mappings);
            haplotype_scorer.decide_barcode_haplotype_support();
            if (haplotype_scorer.barcode_haplotype_mappings.size() > 0) {
                phaseable += 1;
                try {
                    int success = haplotype_scorer.score_haplotypes("formatted_" + output_file);

                    // if we've picked a winner
                    if (success == 0) {
                        std::cout << "Writing output" << std::endl;
                        haplotype_scorer.write_output_success(output_file);
                        graph.write_output_subgraph(haplotype_scorer.winners.first, output_file + ".hap1.fasta", "hap1");
                        graph.write_output_subgraph(haplotype_scorer.winners.second, output_file + ".hap2.fasta", "hap2");
                        phased_success += 1;

                    } else if (success == 1) { // if we're less confident about winner
                        std::cout << "Writing output" << std::endl;

                        haplotype_scorer.write_output_partial_success(output_file);
                        graph.write_output_subgraph(haplotype_scorer.winners.first, "partial_" + output_file + ".hap1.fasta", "hap1");
                        graph.write_output_subgraph(haplotype_scorer.winners.second, "partial_" + output_file + ".hap2.fasta", "hap2");
                        phased_partial_success += 1;

                    }
                } catch (...){
                    std::cout << "Caught exception scoring haplotypes" << std::endl;
                    // no idea what info to put here, apparently can't get full exception from a catch all
                    exceptions += 1;

                }
            }else {
                std::cout << "No mappings suitable for phasing " << std::endl;
                no_mappings += 1;
            }
            std::cout << "----------------------------------------" << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;
        } else {
            unphaseable += 1;
        }
    }
    std::cout << "Phasing " << graphs << " complete, " << phaseable << " contained > 1 bubble, " << unphaseable << " did not." << std::endl;
    std::cout << phased_success << " graphs phased confidently, " << phased_partial_success << " graphs phased less confidently, " << no_mappings << " did not have enough mappings for phasing, and  " << exceptions << " raised." <<std::endl;
    return 0;
}
