//
// Created by Katie Barr (EI) on 15/10/2017.
//

#include "haplotype_scorer.h"


HaplotypeScorer::HaplotypeScorer(std::string mapping_file, std::vector<std::vector <std::string> > possible_haplotypes, Graph graph){
    mapping_filename=mapping_file;
    possible_haplotypes=possible_haplotypes;
    graph = graph;
    for (auto e:graph.edges){
        for (int i=0; i < possible_haplotypes.size(); i++){
                if (std::find(possible_haplotypes[i].begin(), possible_haplotypes[i].end(), e) != possible_haplotypes[i].end()){
                    edge_haplotype_dict[e].push_back(i);
                }
        }
    }
    std::cout << "mappings " <<  mapping_filename <<std::endl;
}

void HaplotypeScorer::add_barcode_vote(std::string barcode, std::string edge, int kmers){
    if (edge_haplotype_dict.find(edge) != edge_haplotype_dict.end()){
        // we only care about mappings to edges in bubbles, which will all have a key in the edge dict
        barcode_edge_mappings[barcode][edge] += kmers;
    }

}
void HaplotypeScorer::load_mappings() {
    std::ifstream infile(mapping_filename);
    std::string line;
    std::string fields[3];
    std::string barcode;
    int counter = 0;
    std::cout << "Loading mappings file " << mapping_filename << std::endl;
    while (std::getline(infile, line)){
        // read name, contig, number of kmers
        std::istringstream(line) >> fields[0] >> fields[1] >> fields[2] ;
        barcode = fields[0].substr(fields[0].find("_") + 1);
        add_barcode_vote(barcode, fields[1], std::stoi(fields[2]));
        counter += 1;

    }
    std::cout << "Loaded " << counter << " mappings" <<std::endl;
}
