//
// Updated by Katie Barr on 11/07/2026.
//

#ifndef PHASER_HAPLOTYPE_SCORER_H
#define PHASER_HAPLOTYPE_SCORER_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <istream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>

#include "graph.h"

class HaplotypeScorer {
private:
    std::vector<std::vector <std::string> > possible_haplotypes;
    std::vector <std::string> unused_barcodes;
    void add_barcode_vote(const std::string&, const std::string&, int);
    // every edge name that appears in any possible haplotype
    std::unordered_set<std::string> haplotype_edges;
    // barcode -> edge name, kmer support
    std::unordered_map<std::string, std::unordered_map<std::string, int> > barcode_edge_mappings;
    std::vector<int> winner_for_barcode(const std::unordered_map<int, int>& haplotype_scores);

public:
    // unordered maps/sets may impact reproducibiilty- if haplotype order mmatters. This was never thecasein our work. 
    std::unordered_map<int, std::unordered_map<std::string, int > > haplotype_barcode_agree;
    std::unordered_map<int, std::unordered_map<std::string, int > > haplotype_barcode_disagree;
    std::unordered_map<std::string, std::unordered_map<int, int> > barcode_haplotype_mappings;
    std::unordered_map<std::string,  int> barcode_hom_mappings;
    std::unordered_map<std::string, int > kmers_per_barcode;
    void write_output_success(std::string);
    void write_output_partial_success(std::string);
    std::unordered_set<std::string> barcodes;
    std::string mapping_filename;
    const Graph& graph;
    HaplotypeScorer(std::string, std::vector<std::vector <std::string> >, Graph&);
    void print_summary(const std::string&, const std::vector<std::pair<int, int> >&, const std::vector<std::pair<int, int> >&, const std::vector<int>&, const std::vector<int>&, const std::vector<int>&);
    //void print_pair_summary(const std::string&, const std::vector<std::pair<std::pair<int, int>, int > >&, const std::vector<std::pair<std::pair<int, int>, int > >&, const std::vector<int>&, const std::vector<int>&, const std::vector<int>& );
    void load_mappings_from_dict(std::map<std::string, std::map<std::string, int> > &);
    void decide_barcode_haplotype_support();
    int max_overall_pair_support;
    int mean_overall_pair_support;
    int max_overall_support;
    int mean_overall_support;
    std::pair<std::vector<std::string>,std::vector<std::string> > winners;
    std::pair<int, int> winning_haplotype;
    int score_haplotypes(std::string);
};


#endif //PHASER_HAPLOTYPE_SCORER_H
