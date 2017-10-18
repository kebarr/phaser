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
    if (bubbles.size() == 0){
        return haplotypes;
    }
    // this assumes all bubbles have 2 contigs, but does not enforce it
    std::vector<std::string> first;
    std::vector<std::string> second;
    first.push_back(std::get<0>(bubbles[0]));
    second.push_back(std::get<1>(bubbles[0]));
    haplotypes.push_back(first);
    haplotypes.push_back(second);
    auto to_index = haplotypes.size();
    auto from_index = haplotypes.size()/2;
    std::string b0;
    std::string b1;
    for (int j=1; j < bubbles.size(); j++){
                b0 = std::get<0>(bubbles[j]);
                to_index = haplotypes.size();
        std::vector<std::vector <std::string> > new_haplotypes;
                for (auto hap: haplotypes) {
                    new_haplotypes.push_back(hap);
                }
                for (auto hap: haplotypes) {
                    new_haplotypes.push_back(hap);
                }
                haplotypes = new_haplotypes;
                for (int i = 0; i < to_index; i++) {
                    haplotypes[i].push_back(b0);
                }
                b1 = std::get<1>(bubbles[j]);

                from_index = (haplotypes.size()) / 2;
                for (int i = from_index; i < haplotypes.size(); i++) {
                    haplotypes[i].push_back(b1);
                }
    }
    return haplotypes;
}


void Graph::load_gfa(std::string infile_name){
    std::ifstream infile(infile_name);
    std::string line;
    std::string fields[5];
    int counter = 0;
    std::cout << "Loading GFA file " << infile_name << std::endl;
    while (std::getline(infile, line)){
        std::istringstream(line) >> fields[0] >> fields[1] >> fields[2] >> fields[3] >> fields[4];
        // to traverse graph only links are required
        if (fields[0] == "L"){
            edges.insert(fields[1]);
            edges.insert(fields[2]);
            std::pair<std::string, std::string> value_fwd = std::make_pair(fields[3], fields[4]);
            // need to store both ways around to ensure every edge connected to a given node is traversed
            std::pair<std::string, std::string> value_bwd = std::make_pair(fields[1], switch_pm[fields[2]]);
            std::pair<std::string, std::string> inverse_link = std::make_pair(fields[3], switch_pm[fields[4]]);
            //std::string t = switch_pm[fields[2]];
            //std::tuple<std::string, std::string> value_bwd = std::make_tuple(t, t);                                                                                                                                                                                                       edge_list[fields[1]].push_back(value_fwd);
            edge_list[std::make_pair(fields[1], fields[2])].insert(value_fwd);
            //std::pair<std::string, std::string> inverse_link = std::make_pair(fields[3], switch_pm[fields[4]]);
            edge_list[inverse_link].insert(value_bwd);
            counter +=1;
        }
    }
    std::cout << "Loaded GFA with " << counter << " links" << edge_list.size()<<std::endl;
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
            if (seq.first != origniating_edge.first){
                return seq;
            }
        }
    }
    return std::make_pair("","");

}

void Graph::traverse_graph(std::string start_node, std::string in_dir, std::vector<std::string > &traversed_edge_list){
    // i replicated the links to ensure every on is a key in the dict- now means we can go same way when supposed to go oppotite ways
    // get nodes joined from other direction- so when we start g
    std::cout << "traversing from " << start_node << std::endl;
    std::pair<std::string, std::string> node = std::make_pair(start_node, switch_pm[in_dir]) ;// should probably just feed this in as aprameter
    std::set<std::pair<std::string, std::string> > adjacent_nodes = edge_list[node];
    std::vector<std::pair<std::string, std::string> > adjacent_nodes_vector;
    std::cout << adjacent_nodes.size() << std::endl;
    for (auto n: adjacent_nodes){
        adjacent_nodes_vector.push_back(n);
    }
    for (auto l:adjacent_nodes){
        std::cout << std::get<0>(l) << " " << std::get<1>(l) << std::endl;

    }
    if (adjacent_nodes.size() == 0){// we can traverse no further
        return;
    } else if (adjacent_nodes.size() == 1 && std::find(traversed_edge_list.begin(), traversed_edge_list.end(), std::get<0>(adjacent_nodes_vector[0])) == traversed_edge_list.end()){
        //travers to next contig
        traversed_edge_list.push_back(std::get<0>(adjacent_nodes_vector[0]));
        traverse_graph(std::get<0>(adjacent_nodes_vector[0]), switch_pm[std::get<1>(adjacent_nodes_vector[0])], traversed_edge_list);
    } else if (std::find(traversed_edge_list.begin(), traversed_edge_list.end(), std::get<0>(adjacent_nodes_vector[0]))== traversed_edge_list.end()
        && std::find(traversed_edge_list.begin(), traversed_edge_list.end(), std::get<0>(adjacent_nodes_vector[1]))== traversed_edge_list.end()){
        traversed_edge_list.push_back(std::get<0>(adjacent_nodes_vector[0]));
        traversed_edge_list.push_back(std::get<0>(adjacent_nodes_vector[1]));

        std::pair<std::string, std::string> contig_other_end_bubble = check_bubble(node, adjacent_nodes_vector);
        if (std::get<0>(contig_other_end_bubble) != "" & std::get<1>(contig_other_end_bubble) != ""){
            //!!!!!! not enforcing bubble degree, but this assumes deg 2....
            //std::cout << "adding bubble " << std::get<0>(adjacent_nodes_vector[0]) << " : " << std::get<0>(adjacent_nodes_vector[1]) << std::endl;
            bubbles.push_back(std::make_pair(std::get<0>(adjacent_nodes_vector[0]), std::get<0>(adjacent_nodes_vector[1])));
                    /// continue traversing from other end of bubble
            traverse_graph(std::get<0>(contig_other_end_bubble), switch_pm[std::get<1>(contig_other_end_bubble)], traversed_edge_list);
        } else {
            // if we've hit something that is not a bubble, we can't phase further, so exit

            return;
        }
    }
}

