//
// Updated by Katie Barr on 11/07/2026.
//

#include <tuple>
#include "graph.h"


Graph::Graph(){
    std::vector<std::vector <int> > bubbles;
}


std::vector<std::vector <std::string> > Graph::calculate_possible_haplotypes(){
    std::vector<std::vector <std::string> > haplotypes;
    if (bubbles.size() == 0){
        return haplotypes;
    }
    // this assumes all bubbles have 2 contigs, but does not enforce it
    std::vector<std::string> first;
    std::vector<std::string> second;
    first.push_back(std::get<0>(bubbles[0])); // start contig of bubble
    second.push_back(std::get<1>(bubbles[0])); // end contig of bubble
    haplotypes.push_back(first);
    haplotypes.push_back(second);
    auto to_index = haplotypes.size(); // number of haplotypes before adding next bubble
    auto from_index = haplotypes.size()/2;
    std::string b0;
    std::string b1;
    for (int j=1; j < bubbles.size(); j++){
                b0 = std::get<0>(bubbles[j]); // name of contig at start of bubble
                to_index = haplotypes.size(); // all the haplotypes pre doubling, add first possible allele to these
        std::vector<std::vector <std::string> > new_haplotypes;
                // duplicate each time eas each possible haplotype must include every possible traversal though the graph
                for (auto hap: haplotypes) {
                    new_haplotypes.push_back(hap);
                }
                for (auto hap: haplotypes) {
                    new_haplotypes.push_back(hap);
                }
                haplotypes = new_haplotypes;
                for (int i = 0; i < to_index; i++) {
                    haplotypes[i].push_back(b0); // start contig for each possible haplotype
                }
                b1 = std::get<1>(bubbles[j]);

                from_index = (haplotypes.size()) / 2; // all the haplotypes post doubling, add second possible allele to these
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
        // GFA is type (link or seq); start contig name; direction; end contig name; direction
        std::istringstream(line) >> fields[0] >> fields[1] >> fields[2] >> fields[3] >> fields[4];
        // to traverse graph only links are required
        if (fields[0] == "L"){
            edges.insert(fields[1]); // start contig
            edges.insert(fields[3]); // end contig
            Strand dir1 = strand_from_gfa(fields[2]);
            Strand dir2 = strand_from_gfa(fields[4]);
            NodeEnd value_fwd{fields[3], dir2};
            // need to store both ways around to ensure every edge connected to a given node is traversed
            NodeEnd value_bwd{fields[1], flip(dir1)};
            NodeEnd inverse_link{fields[3], flip(dir2)};
            edge_list[NodeEnd{fields[1], dir1}].insert(value_fwd);
            edge_list[inverse_link].insert(value_bwd);
            original_edge_dirs[std::make_pair(fields[1], fields[3])] = std::make_pair(dir1, dir2);
            counter +=1;
        } else if (fields[0] == "S"){
            nodes[fields[1]] = fields[2];
        }
    }
    std::cout << "Loaded GFA with " << counter << " links" << edge_list.size()<<std::endl;
}

NodeEnd Graph::check_bubble(NodeEnd origniating_edge, std::vector<NodeEnd> adjacent_nodes){
    // node list are candidate bubble contigs. if the nodes go to and from same contigs, its a bubble
    std::set<NodeEnd> seqs; // use set so same nodes not repeated
    // to be in the same bubble, the contigs have to join the same ends of the adjacent contigs
    for (auto node: adjacent_nodes){  // add each adjecnt node to seqs
        for (auto node2: edge_list[node]) { // check outgoing nodes
            seqs.insert(node2);
        }
        NodeEnd opp_dir_node{node.name, flip(node.strand)};
        for (auto node2: edge_list[opp_dir_node]) { // check incoming nodes
            seqs.insert(node2); // add contig names in other direction
        }
    }
    if (seqs.size() == 2){
        // if only 2 sequences joined to all candidate nodes, they are in a bubble
        // to avoid traversing this part again, return next node and its direction
        for (auto seq: seqs){
            if (seq.name != origniating_edge.name){ //  one of the two elements in seqs is always origniating_edge itself, need this to advance past bubble
                for (auto node: adjacent_nodes){
                    edges_in_bubbles.insert(node.name);
                }
                return seq;
            }
        }
    }
    return NodeEnd{"", Strand::Plus};

}

//TODO: seen at least one example of this stopping one edge earlier than needed
void Graph::traverse_graph(std::string start_node, Strand in_dir, std::vector<std::string > &traversed_edge_list){
    // iterative: traverses graph in specified direction, advancing to the next node when only
    while (true) {
        NodeEnd node{start_node, flip(in_dir)}; /// converts "the end I arrived at start_node through" into "the end I need to leave start_node through," because that's the key edge_list actually indexes on.
        std::set<NodeEnd> adjacent_nodes = edge_list[node]; // all edges collected to this nde
        std::vector<NodeEnd> adjacent_nodes_vector(adjacent_nodes.begin(), adjacent_nodes.end());
        if (adjacent_nodes.size() == 0){// we can traverse no further
            return;
        } else if (adjacent_nodes.size() == 1 && std::find(traversed_edge_list.begin(), traversed_edge_list.end(), adjacent_nodes_vector[0].name) == traversed_edge_list.end()){
            // if the last traversed edge is equal to the first adjacent node
            //travers to next contig
            traversed_edge_list.push_back(adjacent_nodes_vector[0].name);
            start_node = adjacent_nodes_vector[0].name;
            in_dir = flip(adjacent_nodes_vector[0].strand);  // e next contig and the specific end of it that this link attaches to (its arrival end, from the neighbor's perspective). Flipping it before storing it in in_dir keeps the invariant intact for the next lap: at the top of the loop, line 336 will flip it right back, so the next edge_list lookup ends up keyed at exactly the strand value that was just found in
            // the two flips are complementary — one converts "arrival" → "departure" for the current lookup, the other stores next iteration's "arrival" value in a form that, once flipped again, reproduces the departure end correctly.
        } else if (std::find(traversed_edge_list.begin(), traversed_edge_list.end(), adjacent_nodes_vector[0].name)== traversed_edge_list.end()
            && std::find(traversed_edge_list.begin(), traversed_edge_list.end(), adjacent_nodes_vector[1].name)== traversed_edge_list.end()){
            // if there are two adjecent nodes that share the same end point they might be a bubble
            traversed_edge_list.push_back(adjacent_nodes_vector[0].name);
            traversed_edge_list.push_back(adjacent_nodes_vector[1].name);

            NodeEnd contig_other_end_bubble = check_bubble(node, adjacent_nodes_vector);
            if (!contig_other_end_bubble.name.empty()){ // if it is a bubble
                // really lazy but easiest way to check if edge is hom/het for output
                edges_in_bubbles.insert(adjacent_nodes_vector[0].name);
                edges_in_bubbles.insert(adjacent_nodes_vector[1].name);
                bubbles.push_back(std::make_pair(adjacent_nodes_vector[0].name, adjacent_nodes_vector[1].name));
                /// continue traversing from other end of bubble
                start_node = contig_other_end_bubble.name;
                in_dir = flip(contig_other_end_bubble.strand);
            } else {
                // if we've hit something that is not a bubble, we can't phase further, so exit
                return;
            }
        } else {
            // neither branch above matched (e.g. >2 adjacent nodes, or already traversed) - nothing more to do
            return;
        }
    }
}
