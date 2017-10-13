//
// Created by Katie Barr (EI) on 12/10/2017.
//

#include <tuple>
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
    auto to_index = haplotypes.size();
    auto from_index = haplotypes.size();
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

                from_index = (haplotypes.size()) / 2;
                for (int i = 0; i < from_index; i++) {
                    haplotypes[i].push_back(bubble[j]);
                }
            }
        }
    }
    return haplotypes;
}


void Graph::load_gfa(std::string infile_name){
    std::ifstream infile(infile_name);
    std::string line;
    std::string fields[5];
    std::map<std::string, std::string> switch_pm = {{"+","-"}, {"-","+"}};
    while (std::getline(infile, line)){
        std::istringstream(line) >> fields[0] >> fields[1] >> fields[2] >> fields[3] >> fields[4];
        // to traverse graph only links are required
        if (fields[0] == "L"){
            std::tuple<std::string, std::string, std::string> value_fwd = std::make_tuple(fields[3], fields[2], fields[4]);
            // need to store both ways around to ensure every edge connected to a given node is traversed
            std::tuple<std::string, std::string, std::string> value_bwd = std::make_tuple(fields[1], switch_pm[fields[4]], switch_pm[fields[2]]);                                                                                                                                                                                                       edge_list[fields[1]].push_back(value_fwd);
            edge_list[fields[1]].push_back(value_fwd);
            edge_list[fields[3]].push_back(value_bwd);
        }
    }
}

void Graph::traverse_graph(std::string start_node, std::stringstream &link_lines, std::vector<std::pair<std::string, std::string>> &traversed_edge_list){
    std::vector<std::tuple<std::string, std::string, std::string> > adjacent_nodes = edge_list[start_node];
    bool can_traverse_further;// if we can leave an edge in the opposite direction, i.e. if came in on +, and left side of link, then either + on right or - on left
    // now we need to traverse in single direction, so for edge, need to know which direction we came in on  (+/-), and leave by other
    for (auto node = adjacent_nodes.begin(); node != adjacent_nodes.end(); ++node){ //= adjacent_nodes.begin(); node != adjacent_nodes.end(); ++node){
        std::string end_node = std::get<0>(*node);
        std::pair<std::string, std::string> edge_set = std::make_pair(start_node, end_node);
        if (std::find(traversed_edge_list.begin(), traversed_edge_list.end(), edge_set) == traversed_edge_list.end()){
         //determine whether vertex order must be flipped to ensure sequences concatenated in correct order
            bool edge_direction_correct = get<3>(*node);
            if (edge_direction_correct){
                link_lines << "L" << "\t" << start_node << "\t" << std::get<1>(*node)  << "\t" << end_node << "\t" << std::get<2>(*node) << "\t" << "0M" << std::endl;
            } else {
                link_lines << "L" << "\t" << end_node << "\t" << std::get<2>(*node)  << "\t" << start_node << "\t" << std::get<1>(*node) << "\t" << "0M" << std::endl;
            }
            traversed_edge_list.push_back(edge_set);
            if (can_traverse_further){
                traverse_graph(get<0>(*node), link_lines, traversed_edge_list);
            }
    }
}

}
