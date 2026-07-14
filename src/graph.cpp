//
// Updated by Katie Barr on 11/07/2026.
//

#include <tuple>
#include <unordered_set>
#include "graph.h"


Graph::Graph(){
    std::vector<std::vector <int> > bubbles;
}


std::vector<std::vector <std::string> > Graph::calculate_possible_haplotypes(){
    /* Each bubble in a graph increases the number of possible haplotypes by the degree of that bubble
    Because we are only taking actual links from GFAs, every link in a bubble may be part of a haplotype
    Therefore this approach enumerates all of the possible haplotypes arising from graph 
    
    params:  uses attributes on graph already there, no arguments passed
    returns: returns a vector containing each possible haplotype, 
            each haplotype is represented as a vector of strings- the strings being the names of each edge in this haplotype 
            In real applications these strings would be converted to ints to save on memory

    A potential redesign might be to use the barcodes themselves to identify possible haplotypes with minimal support
    */
    std::vector<std::vector <std::string> > haplotypes;
    if (bubbles.size() == 0){
        return haplotypes; // no bubbles, GFA is a single haplotype, no need to enumerate
    }
    haplotypes.push_back(std::vector<std::string>()); // seed with one empty path so the first bubble has something to extend
    int bubble_degree;
    std::string bubble_edge_name;
    for (int j=0; j < bubbles.size(); j++){// for each bubble
        bubble_degree = bubbles[j].size();
        std::vector<std::vector <std::string> > new_haplotypes;
        for (int i = 0; i < bubble_degree; i++){
            bubble_edge_name = bubbles[j][i];
            for (int k = 0; k < haplotypes.size(); k++){
                new_haplotypes.push_back(haplotypes[k]); // replicate existing haplotypes for each possible bubble edge
                new_haplotypes.back().push_back(bubble_edge_name);
            }           
        }
        haplotypes = new_haplotypes;
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
    /*   For each edge that joins to multiple nodes, we need to decide if its a bubble
        A bubble is a set of two edges in a graph that have the same edges leaving the first, 
        and ending at the second. A bubble can have an arbirary number of edges, called the degree.
        In practise, the degreeos of bubbles depend on the ploidy of the organism being assembler 
        and error rate in the sequencing software

        params: originating edge, as type NodeEnc
                adjacent nodes, as a vector of NodeEnd, the edges that join to the originating edge
        returns: The node at the end of the bubble, if the structure is indeed a bubble

        Originally the same set was used for nodes from  both originating edge and adjacent edge 
        This would have silently added non-bubbles as bubbles - the sort of error can equire significant time to notice, 
        possible resulting in incorrect results being presented. For these aplications lab cofirmation is required, 
        it would have been noticed before publication, but this is still a significant losso of lab time. 
        This highlights the need for comprehensive unit testing, particularly for potential edge cases.  
    */
    
    // node list are candidate bubble contigs. if the nodes go to and from same contigs, its a bubble
    std::unordered_set<NodeEnd> seqs_out; // use set so same nodes not repeated, unorderered set because for our applications, haplotype order does not matter. Would need to change this for other use cases. 
    std::unordered_set<NodeEnd> seqs_in;
    // to be in the same bubble, the contigs have to join the same ends of the adjacent contigs
    for (auto node: adjacent_nodes){  // add each adjecnt node to seqs
        auto out_it = edge_list.find(node);
        if (out_it != edge_list.end()) { // check outgoing edges
            for (auto node2: out_it->second) {
                seqs_out.insert(node2);
            }
        }
        NodeEnd opp_dir_node{node.name, flip(node.strand)};
        auto in_it = edge_list.find(opp_dir_node);
        if (in_it != edge_list.end()) { // check incoming edges
            for (auto node2: in_it->second) {
                seqs_in.insert(node2); // add contig names in other direction
            }
        }
    }
    if (seqs_out.size() != 1 || seqs_in.size() != 1){ // if the nodes don't all go to the same contig, its not a bubble- might be part of more complex structure
        return NodeEnd{"", Strand::Plus};
    }
    // to avoid traversing this part again, return next node and its direction
    for (auto seq: seqs_out){
        if (seq.name != origniating_edge.name){ //  one of the two elements in seqs is always origniating_edge itself, need this to advance past bubble
            for (auto node: adjacent_nodes){
                edges_in_bubbles.insert(node.name);
            }
            return seq;
        }
    }
    
    return NodeEnd{"", Strand::Plus};

}

void Graph::traverse_graph(std::string start_node, Strand in_dir){
    /* This function traversees every possible path through the graph defined in the GFA.
       During traversal, it identifies all bubbles, and calculates each possible haplotype
       
         params: start_node, the name of the contig to start traversing from
                in_dir, the direction of the edge we are entering the start_node from

        As nodes can have edges in both directions, this traverses in the spciied direction, the in the opposite
        This ensures all paths are enumerated
    */
    traversed_edge_list.clear(); // each call traverses fresh, independent of any previous call (e.g. the opposite-direction pass)
    // iterative: traverses graph in specified direction, advancing to the next node after determining whether it is a straight link or a bubble
    while (true) {
        NodeEnd node{start_node, flip(in_dir)}; /// converts "the end I arrived at start_node through" into "the end I need to leave start_node through," because that's the key edge_list actually indexes on.
        auto edge_it = edge_list.find(node);
        static const std::set<NodeEnd> empty_node_set;
        const std::set<NodeEnd>& adjacent_nodes = (edge_it != edge_list.end()) ? edge_it->second : empty_node_set; // all edges collected to this node
        std::vector<NodeEnd> adjacent_nodes_vector(adjacent_nodes.begin(), adjacent_nodes.end());
        std::set<std::string> adjacent_node_names;
        for (auto node: adjacent_nodes_vector){
            adjacent_node_names.insert(node.name);
        }
        // check if any ot he nodes has already been traversed
        bool any_already_traversed = std::any_of(adjacent_node_names.begin(), adjacent_node_names.end(), [&](const std::string& n){
            return traversed_edge_list.find(n) != traversed_edge_list.end();
        });
        if (adjacent_nodes.size() == 0){// we can traverse no further
            return;
        } else if (adjacent_nodes.size() == 1 && !any_already_traversed){
            //travers to next contig
            traversed_edge_list.insert(adjacent_nodes_vector[0].name);
            start_node = adjacent_nodes_vector[0].name;
            in_dir = flip(adjacent_nodes_vector[0].strand);  // e next contig and the specific end of it that this link attaches to (its arrival end, from the neighbor's perspective). Flipping it before storing it in in_dir keeps the invariant intact for the next lap: at the top of the loop, line 336 will flip it right back, so the next edge_list lookup ends up keyed at exactly the strand value that was just found in
            // the two flips are complementary — one converts "arrival" → "departure" for the current lookup, the other stores next iteration's "arrival" value in a form that, once flipped again, reproduces the departure end correctly.
        } else if (adjacent_nodes.size() > 1 && !any_already_traversed){
            // if there are several adjacent nodes that share the same end point they might be a (possibly >2-way) bubble
            for (auto& n: adjacent_node_names){
                traversed_edge_list.insert(n);
            }

            NodeEnd contig_other_end_bubble = check_bubble(node, adjacent_nodes_vector);
            if (!contig_other_end_bubble.name.empty()){ // if it is a bubble
                for (auto n: adjacent_nodes_vector){
                    edges_in_bubbles.insert(n.name);
                }
                bubbles.push_back(std::vector<std::string>(adjacent_node_names.begin(), adjacent_node_names.end())); /// add the bubble to the list of bubbles
                /// continue traversing from other end of bubble
                start_node = contig_other_end_bubble.name;
                in_dir = flip(contig_other_end_bubble.strand);
            } else {
                // if we've hit something that is not a bubble, we can't phase further, so exit
                return;
            }
        } else {
            // neither branch above matched (e.g. already traversed) - nothing more to do
            return;
        }
    }
}
