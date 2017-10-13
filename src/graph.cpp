//
// Created by Katie Barr (EI) on 12/10/2017.
//

#include <tuple>
#include "graph.h"


Graph::Graph(){
    std::vector<std::vector <int> > bubbles;
}
/*
void Graph::add_bubble(std::vector <int> bubble){
    bubbles.push_back(bubble);
}*/

std::vector<std::vector <std::string> > Graph::calculate_possible_haplotypes(){
    std::vector<std::vector <std::string> > haplotypes;
    // this assumes all bubbles have 2 contigs, but does not enforce it
    std::vector<std::string> first;
    std::vector<std::string> second;
    first.push_back(std::get<0>(bubbles[0]));
    second.push_back(std::get<1>(bubbles[0]));
    auto to_index = haplotypes.size();
    auto from_index = haplotypes.size();
    std::string b0;
    std::string b1;
    for (auto bubble: bubbles){
                b0 = std::get<0>(bubble);
                to_index = haplotypes.size();
                for (auto &hap: haplotypes) {
                    haplotypes.push_back(hap);
                }
                for (int i = 0; i < to_index; i++) {
                    haplotypes[i].push_back(b0);
                }
                b1 = std::get<1>(bubble);

                from_index = (haplotypes.size()) / 2;
                for (int i = 0; i < from_index; i++) {
                    haplotypes[i].push_back(b1);
                }
    }
    return haplotypes;
}


void Graph::load_gfa(std::string infile_name){
    std::ifstream infile(infile_name);
    std::string line;
    std::string fields[5];
    while (std::getline(infile, line)){
        std::istringstream(line) >> fields[0] >> fields[1] >> fields[2] >> fields[3] >> fields[4];
        // to traverse graph only links are required
        if (fields[0] == "L"){
            std::pair<std::string, std::string> value_fwd = std::make_pair(fields[3], fields[4]);
            // need to store both ways around to ensure every edge connected to a given node is traversed
            std::pair<std::string, std::string> value_bwd = std::make_pair(fields[1], switch_pm[fields[2]]);
            std::pair<std::string, std::string> inverse_link = std::make_pair(fields[3], switch_pm[fields[4]]);
            //std::string t = switch_pm[fields[2]];
            //std::tuple<std::string, std::string> value_bwd = std::make_tuple(t, t);                                                                                                                                                                                                       edge_list[fields[1]].push_back(value_fwd);
            edge_list[std::make_pair(fields[1], fields[2])].push_back(value_fwd);
            //std::pair<std::string, std::string> inverse_link = std::make_pair(fields[3], switch_pm[fields[4]]);
            edge_list[inverse_link].push_back(value_bwd);
        }
    }
}

std::pair<std::string, std::string> Graph::check_bubble(std::pair<std::string, std::string> origniating_edge, std::vector<std::pair<std::string, std::string> > adjacent_nodes){
    // node list are candidate bubble contigs. if the nodes go to and from same contigs, its a bubble
    std::set<std::pair<std::string, std::string> > seqs;
    // to be in the same bubble, the contigs have to join the same ends of the adjacent contigs
    for (auto node: adjacent_nodes){
        for (auto node2: edge_list[node]) {
            seqs.insert(node2);
        }
        std::pair<std::string, std::string>  opp_dir_nodes = std::make_pair(std::get<0>(node), switch_pm[std::get<1>(node)]);
        for (auto node2: edge_list[opp_dir_nodes]) {
            seqs.insert(node2);
        }
    }
    if (seqs.size() == 2){
        // if only 2 sequences joined to all candidate nodes, they are in a bubble
        // to avoid traversing this part again, return next node and its direction
        for (auto seq: seqs){
            if (seq != origniating_edge){
                return seq;
            }
        }
    }
    return std::make_pair("","");

}

void Graph::traverse_graph(std::string start_node, std::string in_dir){
    // get nodes joined from other direction
    std::pair<std::string, std::string> node = std::make_pair(start_node, switch_pm[in_dir]) ;// should probably just feed this in as aprameter
    std::vector<std::pair<std::string, std::string> > adjacent_nodes = edge_list[node];
    if (adjacent_nodes.size() == 0){// we can traverse no further
        exit(0);
    } else if (adjacent_nodes.size() == 1){
        //travers to next contig
    } else {
        std::pair<std::string, std::string> contig_other_end_bubble = check_bubble(node, adjacent_nodes);
        if (std::get<0>(contig_other_end_bubble) != "" & std::get<1>(contig_other_end_bubble) != ""){
            //!!!!!! not enforcing bubble degree, but this assumes deg 2....
            bubbles.push_back(std::make_pair(std::get<0>(adjacent_nodes[0]), std::get<0>(adjacent_nodes[1])));
                    /// continue traversing from other end of bubble
            traverse_graph(std::get<0>(contig_other_end_bubble), std::get<1>(contig_other_end_bubble));
        } else {
            // if we've hit something that is not a bubble, we can't phase further, so exit

            exit(0);
        }
    }
}

