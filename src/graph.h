//
//Updated by Katie Barr on 11/07/2026.
//

#ifndef PHASER_GRAPH_H
#define PHASER_GRAPH_H

#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <istream>
#include <set>
#include <stdlib.h>
#include <algorithm>

enum class Strand {
    Plus,
    Minus
};

inline Strand flip(Strand s) {
    return s == Strand::Plus ? Strand::Minus : Strand::Plus;
}

inline Strand strand_from_gfa(const std::string& s) {
    return s == "+" ? Strand::Plus : Strand::Minus;
}

inline char strand_to_gfa(Strand s) {
    return s == Strand::Plus ? '+' : '-';
}

// a contig plus the strand/end you're entering or leaving it from
struct NodeEnd {
    std::string name;
    Strand strand;

    bool operator==(const NodeEnd& other) const {
        return name == other.name && strand == other.strand;
    }
    bool operator<(const NodeEnd& other) const {
        return name != other.name ? name < other.name : strand < other.strand;
    }
};

namespace std {
    template<>
    struct hash<NodeEnd> {
        size_t operator()(const NodeEnd& n) const noexcept {
            size_t h1 = std::hash<std::string>{}(n.name);
            size_t h2 = std::hash<int>{}(static_cast<int>(n.strand));
            // standard hash-combine: mixes h2 into h1 so name/strand pairs that
            // differ only in strand don't collide in the same bucket
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
        }
    };
}

// it should be possible to extend this so that we can construct phase blocks from haplotypes. keep completely minimal for now.
class Graph
{
private:
    // (contig1, contig2) -> original (dir1, dir2) as read from the GFA link;
    // keyed by name pairs (not NodeEnd) since lookups only have the two names, not their strands
    std::map<std::pair<std::string, std::string>, std::pair<Strand, Strand> > original_edge_dirs;
public:
    void output_contigs_joined_to_contig_list(const std::vector<std::string>&, const std::unordered_map<std::string, int>&, const std::string&) const;
    std::set<std::string> edges;
    std::unordered_map<std::string, std::string>  nodes;
    std::vector<std::vector< std::string> > bubbles;
    std::set<std::string>  edges_in_bubbles;
    NodeEnd check_bubble(NodeEnd, std::vector<NodeEnd> );
    // edge list maps a node end -> set of node ends its connected to
    std::unordered_map <NodeEnd, std::set<NodeEnd> > edge_list;
    void traverse_graph(std::string, Strand, std::set<std::string >&);
    Graph();
    std::vector<std::vector <std::string> >  calculate_possible_haplotypes(void);
    bool can_output_graph_sequence(std::map <NodeEnd, std::vector<NodeEnd> > ) const;
    void load_gfa(std::string);
    // easiest way to actually get phase string is get output sub gfa for each haplotype and stitch together
    void write_output_subgraph(std::vector<std::string> , std::string, std::string  );
    NodeEnd find_start_edge(std::map <NodeEnd, std::vector<NodeEnd> >  ) const;

    void write_sequences_to_file(std::string , std::string,std::vector<std::pair<std::string, bool> > );
};
#endif //PHASER_GRAPH_H
