//
// Created by Katie Barr (EI) on 12/10/2017.
//

#ifndef PHASER_GRAPH_H
#define PHASER_GRAPH_H

#include <vector>

const int BUBBLEDEGREE=2;
// it should be possible to extend this so that we can construct phase blocks from haplotypes. keep completely minimal for now.
class Graph
{
public:
    Graph();
    void add_bubble(std::vector<int>);
    std::vector<std::vector <int> > bubbles;
    std::vector<std::vector <int> > calculate_possible_haplotypes(void);

};
#endif //PHASER_GRAPH_H
#include "graph.cpp"