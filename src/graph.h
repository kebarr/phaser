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

    std::vector<std::pair<std::string, bool> > find_next_edges(std::vector<std::pair<std::string, bool> >, std::vector<std::string> , std::vector<std::string> , std::set<std::pair<std::string, std::string> > );
public:
    void output_contigs_joined_to_contig_list(std::vector<std::string>, std::string);
    std::set<std::string> edges;
    std::map<std::string, std::string>  nodes;
    std::vector<std::pair< std::string, std::string> > bubbles;
    std::set<std::string>  edges_in_bubbles;
    std::pair<std::string, std::string> check_bubble(std::pair<std::string, std::string>, std::vector<std::pair<std::string, std::string> > );
    std::map < std::pair<std::string, std::string> , std::set<std::pair<std::string, std::string> > > edge_list;
    void traverse_graph(std::string, std::string, std::vector<std::string >&);
    Graph();
    std::vector<std::vector <std::string> >  calculate_possible_haplotypes(void);
    bool can_output_graph_sequence(std::map < std::pair<std::string, std::string> , std::vector<std::pair<std::string, std::string> > > );
    void load_gfa(std::string);
    // easiest way to actually get phase string is get output sub gfa for each haplotype and stitch together
    void write_output_subgraph(std::vector<std::string> , std::string, std::string  );
    std::pair<std::string, std::string> find_start_edge(std::map < std::pair<std::string, std::string> , std::vector<std::pair<std::string, std::string> > >  );

    void write_sequences_to_file(std::string , std::string,std::vector<std::pair<std::string, bool> > );
};
#endif //PHASER_GRAPH_H
