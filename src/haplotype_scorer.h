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

#include "graph.h"

class HaplotypeScorer {
private:
    std::vector<std::vector <std::string> > possible_haplotypes;
    void add_barcode_vote(std::string, std::string, int);
    std::map<std::string, std::vector <int> > edge_haplotype_dict;
    std::map<std::string, std::map<std::string, int> > barcode_edge_mappings;

public:
    std::string mapping_filename;
    Graph graph;
    HaplotypeScorer(std::string, std::vector<std::vector <std::string> >, Graph);
    void load_mappings();

};


#endif //PHASER_HAPLOTYPE_SCORER_H