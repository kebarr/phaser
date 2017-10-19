//
// Created by Katie Barr (EI) on 15/10/2017.
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

#include "graph.h"

class HaplotypeScorer {
private:
    std::vector<std::vector <std::string> > possible_haplotypes;
    std::vector <std::string> unused_barcodes;
    void add_barcode_vote(std::string, std::string, int);
    std::map<std::string, std::vector <int> > edge_haplotype_dict;
    std::map<std::string, std::map<std::string, int> > barcode_edge_mappings;
    std::vector<int>  winner_for_barcode(std::string barcode);

public:
    std::map<int, std::map<std::string, int > > haplotype_barcode_agree;
    std::map<int, std::map<std::string, int > > haplotype_barcode_disagree;
    std::map<std::string, std::map<int, int> > barcode_haplotype_mappings;
    std::map<std::string,  int> barcode_hom_mappings;
    std::map<std::string, int > kmers_per_barcode;
    void write_output_success(std::string);
    void write_output_partial_success(std::string);
    std::set <std::string> barcodes;
    std::string mapping_filename;
    Graph graph;
    HaplotypeScorer(std::string, std::vector<std::vector <std::string> >, Graph);
    void load_mappings();

    void load_mappings_from_dict(std::map<std::string, std::map<std::string, int> > &);
    void decide_barcode_haplotype_support();
    int max_overall_pair_support;
    int mean_overall_pair_support;
    int max_overall_support;
    int mean_overall_support;
    std::pair<std::vector<std::string>,std::vector<std::string> > winners;
    std::pair<int, int> winning_pair;
    int score_haplotypes();
};


#endif //PHASER_HAPLOTYPE_SCORER_H
