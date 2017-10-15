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
#include <set>
#include <stdlib.h>
#include <algorithm>


const int BUBBLEDEGREE=2;
// it should be possible to extend this so that we can construct phase blocks from haplotypes. keep completely minimal for now.
class Graph
{
private:
    std::map<std::string, std::string> switch_pm = {{"+","-"}, {"-","+"}};
public:
    std::set<std::string> edges;

    std::vector<std::pair< std::string, std::string> > bubbles;
    std::pair<std::string, std::string> check_bubble(std::pair<std::string, std::string>, std::vector<std::pair<std::string, std::string> > );
    std::map < std::pair<std::string, std::string> , std::set<std::pair<std::string, std::string> > > edge_list;
    void traverse_graph(std::string, std::string, std::vector<std::string >&);
    Graph();
    std::vector<std::vector <std::string> >  calculate_possible_haplotypes(void);
    void load_gfa(std::string);


};
#endif //PHASER_GRAPH_H
