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

std::vector<std::pair<std::string, bool> > Graph::find_next_edges(std::vector<std::pair<std::string, bool> > edges_to_output, std::vector<std::string> edges_seen, std::vector<std::string> bubble_edges, std::set<std::pair<std::string, std::string> > links){
    while(links.size() > 0){
        for (auto link: links) {
            std::string next_edge = std::get<0>(link);
            std::string start_end = std::get<1>(link); // want to go from other end of this seq
            std::string next_dir = switch_pm[start_end];
            links = edge_list[std::make_pair(next_edge, next_dir)];
            if (std::find(edges_seen.begin(), edges_seen.end(), next_edge) == edges_seen.end()) {
                // all links are included in twice to make graph traversal easier - don't want to repeat
                // if edge is in a bubble, it should be in bubble edges
                if (std::find(edges_in_bubbles.begin(), edges_in_bubbles.end(), next_edge) != edges_in_bubbles.end()){
                    if (std::find(bubble_edges.begin(), bubble_edges.end(), next_edge) != bubble_edges.end()){
                        // add sequence to list to output, as we're going from start, it goes at the front
                        // if its a + link, start of next_edge is joined to start of current edge, so reverse it
                        if (start_end == "+"){
                            edges_to_output.insert(edges_to_output.begin(), std::make_pair(next_edge, true));
                        } else {
                            edges_to_output.insert(edges_to_output.begin(), std::make_pair(next_edge, false));
                        }
                    }
                } else {
                    //if its not in a bubble, add it
                    if (start_end == "+"){
                        edges_to_output.insert(edges_to_output.begin(), std::make_pair(next_edge, true));
                    } else {
                        edges_to_output.insert(edges_to_output.begin(), std::make_pair(next_edge, false));
                    }
                }

            }
            edges_seen.push_back(next_edge);
        }
    }
}


void Graph::write_output_subgraph2(std::vector<std::string> bubble_edges, std::string output_file) {
    std::vector<std::string> hom_edges;
    for (auto edge:edges){
        if (std::find(edges_in_bubbles.begin(), edges_in_bubbles.end(), edge) == edges_in_bubbles.end()){
            hom_edges.push_back(edge);
        }
    }
    std::vector < std::pair<std::pair<std::string, std::string> , std::set<std::pair<std::string, std::string> > > > edges_to_include;
    // easier- just go through all links- if its a hom link, or included in bubble edges, take it
    for (auto link:edge_list){
        std::string e1_name = std::get<0>(link.first);
        std::string e2_name = std::get<0>(link.first);
        if (std::find(hom_edges.begin(), hom_edges.end(), e1_name) != hom_edges.end() && std::find(bubble_edges.begin(), bubble_edges.end(), e2_name) != bubble_edges.end()){
            // then this link should be included
            edges_to_include.push_back(link);
        } else if (std::find(hom_edges.begin(), hom_edges.end(), e2_name) != hom_edges.end() && std::find(bubble_edges.begin(), bubble_edges.end(), e1_name) != bubble_edges.end()){
            edges_to_include.push_back(link);

        }
    }
    // to be able to output this as 1 contig, each edge should be joined once at end, once at start - except end ones
    bool can_output = can_output_graph_sequence(edges_to_include);
    std::vector<std::pair<std::string, bool> > edges_to_output;
    if (can_output){
        // need to order/orient contigs - know that apart from ends, each is joined to 1 only at each end
        std::string current_contig =  std::get<0>(std::get<0>(edges_to_include[0]));
        std::string current_contig_dir =  std::get<1>(std::get<0>(edges_to_include[0]));
        edges_to_output.push_back(std::make_pair(current_contig, false));
        // now just follow links from here, in right direction
        std::set<std::pair<std::string, std::string> >  next_contig = std::get<1>(edges_to_include[0]);
        while (next_contig.size() != 0 ){
            for (auto contig: next_contig){
                current_contig = std::get<0>(contig);
                if (std::get<1>(contig) == current_contig_dir){ // then its start to end/end to start- don't reverse
                    edges_to_output.push_back(std::make_pair(current_contig, false));
                } else {
                    edges_to_output.push_back(std::make_pair(current_contig, true));
                }
                // need to go from other end of contig
                current_contig_dir = switch_pm[std::get<1>(contig)];

            }
            next_contig.clear();
            for (auto l: edges_to_include){
                if (std::get<0>(l) == std::make_pair(current_contig, current_contig_dir)){
                    next_contig= std::get<1>(l);
                }
            }
        }
        std::string opposite_dir = switch_pm[std::get<1>(std::get<0>(edges_to_include[0]))];
        next_contig.clear();
        for (auto l: edges_to_include){
            if (std::get<0>(l) == std::make_pair(std::get<0>(std::get<0>(edges_to_include[0])), opposite_dir)){
                next_contig= std::get<1>(l);
            }
        }       

    }
}

bool Graph::can_output_graph_sequence(std::vector < std::pair<std::pair<std::string, std::string> , std::set<std::pair<std::string, std::string> > > >  edges){
    std::map<std::string, std::set<std::string> > edges_start;
    std::map<std::string, std::set<std::string> > edges_end;
    // edge dict replicates links- have from_link, from_start_end : to_link to_start_end
    // avoid repetition by only going through dict keys
    // need each edge joined to one contig at start, one contig at end
    for (auto link:edges){
        for (auto linked_to: std::get<1>(link)) {
            if (std::get<1>(std::get<0>(link)) == "+") {
                edges_end[std::get<0>(std::get<0>(link))].insert(std::get<0>(linked_to));
            }
            if (std::get<1>(std::get<0>(link)) == "-") {
                edges_start[std::get<0>(std::get<0>(link))].insert(std::get<0>(linked_to));
            }
        }

    }
    for (auto e: edges_end){
        if (e.second.size() > 1){
            return false;
        }
    }
    for (auto e: edges_start){
        if (e.second.size() > 1){
            return false;
        }
    }
    return true;
}


void Graph::write_output_subgraph(std::vector<std::string> bubble_edges, std::string output_file) {
    // bubble edges are the edges we want to include, in weird order due to traversing from middle outwards
    // for each bubble edge, find edge its attached to front/back and build up sequence
    // as bubbles aren't guaranteed adjacent, will build up parts then stitch together
    std::vector<std::string> edges_seen;
    // these should be in orde, bool is whether sequence is reversed or not
    std::vector<std::pair<std::string, bool> > edges_to_output;
    std::string b = bubble_edges[0];
    // build up list of edges to link in end - start order -
    // start at arbitrary bubble, traverse from there building up list in order
    edges_seen.push_back(b);
    // assume we have beautiful bubble-contig-bubble-structure for now
    // so these are both length 0 or 1
    std::set<std::pair<std::string, std::string> > links_start = edge_list[std::make_pair(b, "-")];
    std::set<std::pair<std::string, std::string> > links_end = edge_list[std::make_pair(b, "+")];
    std::string seq = nodes[b];
    // should be max 1 each of links start and end, if we've reached end, can't go further



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
        } else if (fields[0] == "S"){
            nodes[fields[2]] = fields[3];
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
            // really lazy but easiest way to check if edge is hom/het for output
            edges_in_bubbles.push_back(std::get<0>(adjacent_nodes_vector[0]));
            edges_in_bubbles.push_back(std::get<0>(adjacent_nodes_vector[1]));
            bubbles.push_back(std::make_pair(std::get<0>(adjacent_nodes_vector[0]), std::get<0>(adjacent_nodes_vector[1])));
                    /// continue traversing from other end of bubble
            traverse_graph(std::get<0>(contig_other_end_bubble), switch_pm[std::get<1>(contig_other_end_bubble)], traversed_edge_list);
        } else {
            // if we've hit something that is not a bubble, we can't phase further, so exit

            return;
        }
    }
}

