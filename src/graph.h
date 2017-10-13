//
// Created by Katie Barr (EI) on 12/10/2017.
//

#ifndef PHASER_GRAPH_H
#define PHASER_GRAPH_H

#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <istream>


const int BUBBLEDEGREE=2;
// it should be possible to extend this so that we can construct phase blocks from haplotypes. keep completely minimal for now.
class Graph
{
public:
    std::map < std::pair<std::string, std::string> , std::vector<std::pair<std::string, std::string> > > edge_list;
    void traverse_graph(std::string, std::stringstream&, std::vector<std::pair<std::string, std::string>>&, std::string);
    Graph();
    void add_bubble(std::vector<int>);
    std::vector<std::vector <int> > bubbles;
    std::vector<std::vector <int> > calculate_possible_haplotypes(void);
    void load_gfa(std::string);


};
#endif //PHASER_GRAPH_H
