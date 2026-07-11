//
//Updated by Katie Barr on 11/07/2026.
//
// output/serialization logic split out of graph.cpp: reconstructing a linear
// contig path (for a chosen haplotype or a homozygous subgraph) and writing it out.

#include <tuple>
#include <unordered_set>
#include "graph.h"


void Graph::write_output_subgraph(std::vector<std::string> bubble_edges, std::string output_file, std::string sequence_name) {
    std::unordered_set<std::string> hom_edges;
    for (auto edge:edges){
        if (edges_in_bubbles.find(edge) == edges_in_bubbles.end()){
            hom_edges.insert(edge);
        }
    }
    std::unordered_set<std::string> bubble_edge_set(bubble_edges.begin(), bubble_edges.end());
    std::map <NodeEnd, std::vector<NodeEnd> > edges_to_include;
    // easier- just go through all links- if its a hom link, or included in bubble edges, take it
    for (auto link:edge_list){
        std::string e1_name = link.first.name;
        for (auto joined_to: link.second) {
            std::string e2_name = joined_to.name;
            if (hom_edges.count(e1_name) && bubble_edge_set.count(e2_name)) {
                // then this link should be included
                edges_to_include[link.first].push_back(joined_to);
            } else if (hom_edges.count(e2_name) && bubble_edge_set.count(e1_name)) {
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

void Graph::output_contigs_joined_to_contig_list(std::vector<std::string> bubble_edges, std::map<std::string, int > agreeing_barcodes, std::string outfile_name) const{
    // need to be able to reconstruct each haplotype sequence with its winning contig, phaser just outputs contig choices, so need inbetween links
    std::unordered_set<std::string> hom_edges;
    for (auto edge:edges){
        if (edges_in_bubbles.find(edge) == edges_in_bubbles.end()){
            hom_edges.insert(edge); // all homozygous edges
        }
    }
    std::unordered_set<std::string> bubble_edge_set(bubble_edges.begin(), bubble_edges.end());
    std::map <NodeEnd, std::vector<NodeEnd> > edges_to_include;
    // easier- just go through all links- if its a hom link, or included in bubble edges, take it
    for (auto link:edge_list){
        std::string e1_name = link.first.name;
        for (auto joined_to: link.second) {
            std::string e2_name = joined_to.name;
            // keep a link only if it connects a homozygous contig to this haplotype's chosen bubble allele — i.e. the boundary edges where a hom stretch meets the specific allele this haplotype picked. That's needed because possible_haplotypes only lists the bubble alleles in isolation; to actually reconstruct a full contig path you need to find what hom contig each chosen allele attaches to on either side.
            if (hom_edges.count(e1_name) && bubble_edge_set.count(e2_name)) {
                // if first contig is not in homologous edges and second contig is not in bubble edges join them
                // then this link should be included
                edges_to_include[link.first].push_back(joined_to);
            } else if (hom_edges.count(e2_name) && bubble_edge_set.count(e1_name)) {
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
        while (edge_leaving_other_way.size() != 0) { // traverse until no more edges in this direction
            edges_to_output.push_back(edge_name);
            edge_leaving_other_way = edges_to_include[NodeEnd{edge_name, current_dir}];
            next_edge = edge_leaving_other_way;
            edge_name = next_edge[0].name;
            current_dir = next_edge[0].strand;
        }
        std::ofstream out(outfile_name);
        std::cout << "Writing output to  " << outfile_name <<std::endl;;
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
        out  <<std::endl;;
        for (auto b:agreeing_barcodes){
            out << b.first  << std::endl;;
        }

    }
}
