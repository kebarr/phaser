//
// Created by Katie Barr (EI) on 12/10/2017.
//

#include "graph.h"


Graph::Graph(){
    std::vector<std::vector <int> > bubbles;
}

void Graph::add_bubble(std::vector <int> bubble){
    bubbles.push_back(bubble);
}

std::vector<std::vector <int> > Graph::calculate_possible_haplotypes(){
    std::vector<std::vector <int> > haplotypes;
    // this assumes all bubbles have 2 contigs, but does not enforce it
    std::vector<int> first;
    std::vector<int> second;
    first.push_back(bubbles[0][0]);
    second.push_back(bubbles[0][1]);
    auto to_index;
    auto from_index;
    for (auto bubble: bubbles){
        for (int j=0; bubble.size(); j++) {
            if (j == 0) {
                to_index = haplotypes.size();
                for (auto &hap: haplotypes) {
                    haplotypes.push_back(hap);
                }
                for (int i = 0; i < to_index; i++) {
                    haplotypes[i].push_back(bubble[j]);
                }
            } else {

                from_index = (haplotypes.size()) / 2
                for (int i = 0; i < from_index; i++) {
                    haplotypes[i].push_back(bubble[j]);
                }
            }
        }
    }
    return haplotypes;
}