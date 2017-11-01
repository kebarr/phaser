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


void Graph::write_output_subgraph(std::vector<std::string> bubble_edges, std::string output_file, std::string sequence_name) {
    std::vector<std::string> hom_edges;
    for (auto edge:edges){
        if (std::find(edges_in_bubbles.begin(), edges_in_bubbles.end(), edge) == edges_in_bubbles.end()){
            hom_edges.push_back(edge);
        }
    }
    std::map < std::pair<std::string, std::string> , std::vector<std::pair<std::string, std::string> > > edges_to_include;
    // easier- just go through all links- if its a hom link, or included in bubble edges, take it
    for (auto link:edge_list){
        std::string e1_name = std::get<0>(link.first);
        for (auto joined_to: link.second) {
            std::string e2_name = std::get<0>(joined_to);
            if (std::find(hom_edges.begin(), hom_edges.end(), e1_name) != hom_edges.end() &&
                std::find(bubble_edges.begin(), bubble_edges.end(), e2_name) != bubble_edges.end()) {
                // then this link should be included
                edges_to_include[link.first].push_back(joined_to);
            } else if (std::find(hom_edges.begin(), hom_edges.end(), e2_name) != hom_edges.end() &&
                       std::find(bubble_edges.begin(), bubble_edges.end(), e1_name) != bubble_edges.end()) {
                edges_to_include[link.first].push_back(joined_to);

            }
        }
    }
    // to be able to output this as 1 contig, each edge should be joined once at end, once at start - except end ones
    bool can_output = can_output_graph_sequence(edges_to_include);
    std::vector<std::pair<std::string, bool> > edges_to_output;
    if (can_output){
        // need to order/orient contigs - know that apart from ends, each is joined to 1 only at each end
        // ok, try again, find one of end contigs and just go along
        auto start_edge = find_start_edge(edges_to_include);
        auto previous_dir = std::get<1>(start_edge);
        std::vector<std::pair<std::string, std::string> > next_edge = edges_to_include[std::make_pair(std::get<0>(start_edge), previous_dir)];
        edges_to_output.push_back(std::make_pair(std::get<0>(start_edge), false));
        auto e =  next_edge[0];
        auto edge_name = std::get<0>(e);
        std::string current_dir = std::get<1>(e);
        auto edge_leaving_other_way = edges_to_include[std::make_pair(edge_name, switch_pm[current_dir])];
        while (edge_leaving_other_way.size() != 0) {
            if (current_dir == previous_dir) {
                edges_to_output.push_back(std::make_pair(edge_name, false));
            } else {
                    edges_to_output.push_back(std::make_pair(edge_name, true));

            }
            edge_leaving_other_way = edges_to_include[std::make_pair(edge_name, current_dir)];
            next_edge = edge_leaving_other_way;
            edge_name = std::get<0>(next_edge[0]);
            previous_dir = current_dir;
            current_dir = std::get<1>(next_edge[0]);
            }
        write_sequences_to_file(output_file, sequence_name, edges_to_output);


    }
}

void Graph::write_sequences_to_file(std::string output_filename, std::string sequence_name, std::vector<std::pair<std::string, bool> > edges_to_output){
    std::string sequence;
    for (auto edge: edges_to_output){
        auto seq = nodes[std::get<0>(edge)];
        if (std::get<1>(edge)){
            std::reverse(seq.begin(), seq.end());
        }
        sequence = sequence + seq;
    }
    std::ofstream out(output_filename);
    out << ">" << sequence_name << std::endl << sequence << std::endl;
}

std::pair<std::string, std::string> Graph::find_start_edge(std::map < std::pair<std::string, std::string> , std::vector<std::pair<std::string, std::string> > >  edges_to_subgraph){
    for (auto e: edges_to_subgraph){
        std::pair<std::string, std::string> inverse_links = std::make_pair(std::get<0>(e.first), switch_pm[std::get<1>(e.first)]);
        if (e.second.size() == 0 or edges_to_subgraph.find(inverse_links) == edges_to_subgraph.end()){
            std::pair<std::string, std::string> start_edge = std::make_pair(std::get<0>(e.first), std::get<1>(e.first));
            return start_edge;
        }
    }
    return std::make_pair("","");
}

bool Graph::can_output_graph_sequence(std::map < std::pair<std::string, std::string> , std::vector<std::pair<std::string, std::string> > >  edges){
    std::map<std::string, std::set<std::string> > edges_start;
    std::map<std::string, std::set<std::string> > edges_end;
    // edge dict replicates links- have from_link, from_start_end : to_link to_start_end
    // avoid repetition by only going through dict keys
    // need each edge joined to one contig at start, one contig at end
    for (auto link:edges){
        for (auto linked_to: link.second) {
            if (std::get<1>(link.first) == "+") {// links joined to the end of this
                std::string edge_from = std::get<0>(link.first);
                edges_end[edge_from].insert(std::get<0>(linked_to));
            }
            if (std::get<1>(link.first) == "-") { // links joined to start of this - so go before it in list
                std::string edge_from = std::get<0>(link.first);
                edges_start[edge_from].insert(std::get<0>(linked_to));

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
            edges.insert(fields[3]);
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
            nodes[fields[1]] = fields[2];
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
                for (auto node: adjacent_nodes){
                    edges_in_bubbles.insert(std::get<0>(node));
                }
                return seq;
            }
        }
    }
    return std::make_pair("","");

}

void Graph::output_contigs_joined_to_contig_list(std::vector<std::string> bubble_edges, std::string outfile_name){
    // need to be able to reconstruct each haplotype sequence, phaser just outputs contig chcoices, so need inbetween links
    std::vector<std::string> hom_edges;
    for (auto edge:edges){
        if (std::find(edges_in_bubbles.begin(), edges_in_bubbles.end(), edge) == edges_in_bubbles.end()){
            hom_edges.push_back(edge);
        }
    }
    std::map < std::pair<std::string, std::string> , std::vector<std::pair<std::string, std::string> > > edges_to_include;
    // easier- just go through all links- if its a hom link, or included in bubble edges, take it
    for (auto link:edge_list){
        std::string e1_name = std::get<0>(link.first);
        for (auto joined_to: link.second) {
            std::string e2_name = std::get<0>(joined_to);
            if (std::find(hom_edges.begin(), hom_edges.end(), e1_name) != hom_edges.end() &&
                std::find(bubble_edges.begin(), bubble_edges.end(), e2_name) != bubble_edges.end()) {
                // then this link should be included
                edges_to_include[link.first].push_back(joined_to);
            } else if (std::find(hom_edges.begin(), hom_edges.end(), e2_name) != hom_edges.end() &&
                       std::find(bubble_edges.begin(), bubble_edges.end(), e1_name) != bubble_edges.end()) {
                edges_to_include[link.first].push_back(joined_to);

            }
        }
    }
    // to be able to output this as 1 contig, each edge should be joined once at end, once at start - except end ones
    bool can_output = can_output_graph_sequence(edges_to_include);
    std::vector<std::string > edges_to_output;
    if (can_output) {
        // need to order/orient contigs - know that apart from ends, each is joined to 1 only at each end
        // ok, try again, find one of end contigs and just go along
        auto start_edge = find_start_edge(edges_to_include);
        auto previous_dir = std::get<1>(start_edge);
        std::vector<std::pair<std::string, std::string> > next_edge = edges_to_include[std::make_pair(
                std::get<0>(start_edge), previous_dir)];
        edges_to_output.push_back(std::get<0>(start_edge));
        auto e = next_edge[0];
        auto edge_name = std::get<0>(e);
        std::string current_dir = std::get<1>(e);
        auto edge_leaving_other_way = edges_to_include[std::make_pair(edge_name, switch_pm[current_dir])];
        while (edge_leaving_other_way.size() != 0) {
            edges_to_output.push_back(edge_name);
            edge_leaving_other_way = edges_to_include[std::make_pair(edge_name, current_dir)];
            next_edge = edge_leaving_other_way;
            edge_name = std::get<0>(next_edge[0]);
            current_dir = std::get<1>(next_edge[0]);
        }
        std::ofstream out(outfile_name);
        for (auto edge:edges_to_output){
            out << edge << "\t";
        }
        out << "\n";

    }
}

void Graph::traverse_graph(std::string start_node, std::string in_dir, std::vector<std::string > &traversed_edge_list){
    // i replicated the links to ensure every on is a key in the dict- now means we can go same way when supposed to go oppotite ways
    // get nodes joined from other direction- so when we start g
    std::pair<std::string, std::string> node = std::make_pair(start_node, switch_pm[in_dir]) ;// should probably just feed this in as aprameter
    std::set<std::pair<std::string, std::string> > adjacent_nodes = edge_list[node];
    std::vector<std::pair<std::string, std::string> > adjacent_nodes_vector;
    for (auto n: adjacent_nodes){
        adjacent_nodes_vector.push_back(n);
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
            edges_in_bubbles.insert(std::get<0>(adjacent_nodes_vector[0]));
            edges_in_bubbles.insert(std::get<0>(adjacent_nodes_vector[1]));
            bubbles.push_back(std::make_pair(std::get<0>(adjacent_nodes_vector[0]), std::get<0>(adjacent_nodes_vector[1])));
                    /// continue traversing from other end of bubble
            traverse_graph(std::get<0>(contig_other_end_bubble), switch_pm[std::get<1>(contig_other_end_bubble)], traversed_edge_list);
        } else {
            // if we've hit something that is not a bubble, we can't phase further, so exit

            return;
        }
    }
}

