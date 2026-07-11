//
// Created by Katie Barr (EI) on 12/10/2017.
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

std::vector<std::pair<std::string, bool> > Graph::find_next_edges(std::vector<std::pair<std::string, bool> > edges_to_output, std::vector<std::string> edges_seen, std::vector<std::string> bubble_edges, std::set<NodeEnd> links){
    while(links.size() > 0){
        std::set<NodeEnd> next_links;
        for (auto link: links) {
            std::string next_edge = link.name;
            Strand start_end = link.strand; // want to go from other end of this seq
            next_links = edge_list[NodeEnd{next_edge, flip(start_end)}];
            if (std::find(edges_seen.begin(), edges_seen.end(), next_edge) == edges_seen.end()) {
                // all links are included in twice to make graph traversal easier - don't want to repeat
                // if edge is in a bubble, it should be in bubble edges
                if (std::find(edges_in_bubbles.begin(), edges_in_bubbles.end(), next_edge) != edges_in_bubbles.end()){
                    if (std::find(bubble_edges.begin(), bubble_edges.end(), next_edge) != bubble_edges.end()){
                        // add sequence to list to output, as we're going from start, it goes at the front
                        // if its a + link, start of next_edge is joined to start of current edge, so reverse it
                        if (start_end == Strand::Plus){
                            edges_to_output.insert(edges_to_output.begin(), std::make_pair(next_edge, true));
                        } else {
                            edges_to_output.insert(edges_to_output.begin(), std::make_pair(next_edge, false));
                        }
                    }
                } else {
                    //if its not in a bubble, add it
                    if (start_end == Strand::Plus){
                        edges_to_output.insert(edges_to_output.begin(), std::make_pair(next_edge, true));
                    } else {
                        edges_to_output.insert(edges_to_output.begin(), std::make_pair(next_edge, false));
                    }
                }

            }
            edges_seen.push_back(next_edge);
        }
        links = next_links;
    }
    return edges_to_output;
}


void Graph::write_output_subgraph(std::vector<std::string> bubble_edges, std::string output_file, std::string sequence_name) {
    std::vector<std::string> hom_edges;
    for (auto edge:edges){
        if (std::find(edges_in_bubbles.begin(), edges_in_bubbles.end(), edge) == edges_in_bubbles.end()){
            hom_edges.push_back(edge);
        }
    }
    std::map <NodeEnd, std::vector<NodeEnd> > edges_to_include;
    // easier- just go through all links- if its a hom link, or included in bubble edges, take it
    for (auto link:edge_list){
        std::string e1_name = link.first.name;
        for (auto joined_to: link.second) {
            std::string e2_name = joined_to.name;
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
        auto previous_dir = start_edge.strand;
        std::vector<NodeEnd> next_edge = edges_to_include[NodeEnd{start_edge.name, previous_dir}];
        edges_to_output.push_back(std::make_pair(start_edge.name, false));
        auto e =  next_edge[0];
        auto edge_name = e.name;
        Strand current_dir = e.strand;
        auto edge_leaving_other_way = edges_to_include[NodeEnd{edge_name, flip(current_dir)}];
        while (edge_leaving_other_way.size() != 0) {
            if (current_dir == previous_dir) {
                edges_to_output.push_back(std::make_pair(edge_name, false));
            } else {
                    edges_to_output.push_back(std::make_pair(edge_name, true));

            }
            edge_leaving_other_way = edges_to_include[NodeEnd{edge_name, current_dir}];
            next_edge = edge_leaving_other_way;
            edge_name = next_edge[0].name;
            previous_dir = current_dir;
            current_dir = next_edge[0].strand;
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

NodeEnd Graph::find_start_edge(std::map <NodeEnd, std::vector<NodeEnd> >  edges_to_subgraph) const{
    for (auto e: edges_to_subgraph){
        NodeEnd inverse_link{e.first.name, flip(e.first.strand)};
        if (e.second.size() == 0 or edges_to_subgraph.find(inverse_link) == edges_to_subgraph.end()){
            return e.first;
        }
    }
    return NodeEnd{"", Strand::Plus};
}

bool Graph::can_output_graph_sequence(std::map <NodeEnd, std::vector<NodeEnd> >  edges) const{
    std::map<std::string, std::set<std::string> > edges_start;
    std::map<std::string, std::set<std::string> > edges_end;
    // edge dict replicates links- have from_link, from_start_end : to_link to_start_end
    // avoid repetition by only going through dict keys
    // need each edge joined to one contig at start, one contig at end
    for (auto link:edges){
        for (auto linked_to: link.second) {
            if (link.first.strand == Strand::Plus) {// links joined to the end of this
                edges_end[link.first.name].insert(linked_to.name);
            } else { // links joined to start of this - so go before it in list
                edges_start[link.first.name].insert(linked_to.name);
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
    std::set<NodeEnd> seqs;
    // to be in the same bubble, the contigs have to join the same ends of the adjacent contigs
    for (auto node: adjacent_nodes){
        for (auto node2: edge_list[node]) {
            seqs.insert(node2);
        }
        NodeEnd opp_dir_node{node.name, flip(node.strand)};
        for (auto node2: edge_list[opp_dir_node]) {
            seqs.insert(node2);
        }
    }
    if (seqs.size() == 2){
        // if only 2 sequences joined to all candidate nodes, they are in a bubble
        // to avoid traversing this part again, return next node and its direction
        for (auto seq: seqs){
            if (seq.name != origniating_edge.name){
                for (auto node: adjacent_nodes){
                    edges_in_bubbles.insert(node.name);
                }
                return seq;
            }
        }
    }
    return NodeEnd{"", Strand::Plus};

}

void Graph::output_contigs_joined_to_contig_list(std::vector<std::string> bubble_edges, std::map<std::string, int > agreeing_barcodes, std::string outfile_name) const{
    // need to be able to reconstruct each haplotype sequence, phaser just outputs contig chcoices, so need inbetween links
    std::vector<std::string> hom_edges;
    for (auto edge:edges){
        if (std::find(edges_in_bubbles.begin(), edges_in_bubbles.end(), edge) == edges_in_bubbles.end()){
            hom_edges.push_back(edge);
        }
    }
    std::map <NodeEnd, std::vector<NodeEnd> > edges_to_include;
    // easier- just go through all links- if its a hom link, or included in bubble edges, take it
    for (auto link:edge_list){
        std::string e1_name = link.first.name;
        for (auto joined_to: link.second) {
            std::string e2_name = joined_to.name;
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
        auto previous_dir = start_edge.strand;
        std::vector<NodeEnd> next_edge = edges_to_include[NodeEnd{start_edge.name, previous_dir}];
        edges_to_output.push_back(start_edge.name);
        auto e = next_edge[0];
        auto edge_name = e.name;
        Strand current_dir = e.strand;
        auto edge_leaving_other_way = edges_to_include[NodeEnd{edge_name, flip(current_dir)}];
        while (edge_leaving_other_way.size() != 0) {
            edges_to_output.push_back(edge_name);
            edge_leaving_other_way = edges_to_include[NodeEnd{edge_name, current_dir}];
            next_edge = edge_leaving_other_way;
            edge_name = next_edge[0].name;
            current_dir = next_edge[0].strand;
        }
        std::ofstream out(outfile_name);
        for (auto i=0;i < edges_to_output.size() -1; i++){
            auto edge = edges_to_output[i];
            auto next_edge = edges_to_output[i+1];
            if (original_edge_dirs.find(std::make_pair(edge, next_edge)) != original_edge_dirs.end()){
                Strand dir = original_edge_dirs.at(std::make_pair(edge, next_edge)).first;
                out << edge << strand_to_gfa(dir) << ",";
            } else if (original_edge_dirs.find(std::make_pair(next_edge, edge)) != original_edge_dirs.end()){
                // if original link was in opposite dir, need to switch plus/minus
                Strand dir = flip(original_edge_dirs.at(std::make_pair(next_edge, edge)).first);
                out << edge << strand_to_gfa(dir) << ",";
            }
        }
        out << "\n";
        for (auto b:agreeing_barcodes){
            out << b.first << "\n";
        }

    }
}

//TODO: seen at least one example of this stopping one edge earlier than needed
void Graph::traverse_graph(std::string start_node, Strand in_dir, std::vector<std::string > &traversed_edge_list){
    // iterative: traverses graph in specified direction, advancing to the next node when only
    // 1 node is found (so no phasing required). written as a loop rather than recursion since
    // every recursive call here was a tail call and a long contig chain could otherwise overflow the stack
    // i replicated the links to ensure every on is a key in the dict- now means we can go same way when supposed to go oppotite ways
    // get nodes joined from other direction- so when we start g
    while (true) {
        NodeEnd node{start_node, flip(in_dir)};
        std::set<NodeEnd> adjacent_nodes = edge_list[node]; // all edges collected to this nde
        std::vector<NodeEnd> adjacent_nodes_vector(adjacent_nodes.begin(), adjacent_nodes.end());
        if (adjacent_nodes.size() == 0){// we can traverse no further
            return;
        } else if (adjacent_nodes.size() == 1 && std::find(traversed_edge_list.begin(), traversed_edge_list.end(), adjacent_nodes_vector[0].name) == traversed_edge_list.end()){
            //travers to next contig
            traversed_edge_list.push_back(adjacent_nodes_vector[0].name);
            start_node = adjacent_nodes_vector[0].name;
            in_dir = flip(adjacent_nodes_vector[0].strand);
        } else if (std::find(traversed_edge_list.begin(), traversed_edge_list.end(), adjacent_nodes_vector[0].name)== traversed_edge_list.end()
            && std::find(traversed_edge_list.begin(), traversed_edge_list.end(), adjacent_nodes_vector[1].name)== traversed_edge_list.end()){
            // if there are two adjecent nodes that share the same end point they might be a bubble
            traversed_edge_list.push_back(adjacent_nodes_vector[0].name);
            traversed_edge_list.push_back(adjacent_nodes_vector[1].name);

            NodeEnd contig_other_end_bubble = check_bubble(node, adjacent_nodes_vector);
            if (!contig_other_end_bubble.name.empty()){
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
