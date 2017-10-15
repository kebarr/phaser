//
// Created by Katie Barr (EI) on 15/10/2017.
//

#include "haplotype_scorer.h"


HaplotypeScorer::HaplotypeScorer(std::string mapping_file, std::vector<std::vector <std::string> > possible_hs, Graph g){
    mapping_filename=mapping_file;
    possible_haplotypes=possible_hs;
    graph = g;
    for (auto e:graph.edges){
        for (int i=0; i < possible_haplotypes.size(); i++){
                if (std::find(possible_haplotypes[i].begin(), possible_haplotypes[i].end(), e) != possible_haplotypes[i].end()){
                    edge_haplotype_dict[e].push_back(i);
                }
        }
    }
    std::cout << "mappings " <<  mapping_filename <<std::endl;
    std::cout << "edges in haplotypes: " << edge_haplotype_dict.size() <<std::endl;
    for (auto e:edge_haplotype_dict){
        std::cout << e.first << " ";
    }
    std::cout << std::endl;
}

void print_vector(std::vector<std::string> vec){
    for (auto a: vec){
        std::cout << a << " " << std::endl;
    }
    std::cout << std::endl;
}

void HaplotypeScorer::decide_barcode_haplotype_support(){
    std::vector<std::string> edges;
    std::vector<std::string> edges_in_haplotype;
    std::vector<std::string> h;
    int support;
    int haplotypes_supported = 0;
    std::cout << "Calculating barcode haplotype support" << std::endl;
    for (auto mapping:barcode_edge_mappings){
        std::cout << "Checking barcode " << mapping.first <<std::endl;
        // if barcode maps to more than 1 edge in bubbles
        if (mapping.second.size() > 1){
            for (auto e: mapping.second){
                edges.push_back(e.first);
            }
            //std::cout << "edges barcode maps to: " << edges.size();
            //print_vector(edges);
            for (int i=0; i < possible_haplotypes.size(); i++){
                h = possible_haplotypes[i];
                //std::cout << "edges in haplotype: " << h.size();
                //print_vector(h);
                // find all edges in each haplotype that this barcode maps to
                for (auto e1: h){
                    //std::cout << e1 << std::endl;
                    if (std::find(edges.begin(), edges.end(), e1) != edges.end()){
                        edges_in_haplotype.push_back(e1);
                        //std::cout << "added "<< e1 << std::endl;
                    }
                }

                //std::cout << "edges in haplotype barcode maps to: " << edges.size();
                //print_vector(edges_in_haplotype);
                // somewhat arbitrary rule to decide if the barcode supports a haplotype enough
                if (edges_in_haplotype.size() > edges.size()/2 && edges_in_haplotype.size() >1){
                    for (auto a: edges_in_haplotype){
                        support += mapping.second[a];
                    }
                    barcode_haplotype_mappings[mapping.first][i] = support;
                    support = 0;
                    haplotypes_supported +=1;
                }
                edges_in_haplotype.clear();
            }

            edges.clear();

        } else {
            unused_barcodes.push_back(mapping.first);
        }
        std::cout << "barcode " << mapping.first << " supports " << haplotypes_supported << std::endl;

        haplotypes_supported = 0;
    }
    std::cout << "Calculated haplotype support for each barcode, " << barcode_haplotype_mappings.size() << " useful barcodes " <<std::endl;

}

void HaplotypeScorer::add_barcode_vote(std::string barcode, std::string edge, int kmers){
    if (edge_haplotype_dict.find(edge) != edge_haplotype_dict.end()){
        // we only care about mappings to edges in bubbles, which will all have a key in the edge dict
        barcode_edge_mappings[barcode][edge] += kmers;
    }

}
void HaplotypeScorer::load_mappings() {
    std::ifstream infile(mapping_filename);
    std::string line;
    std::string fields[3];
    std::string barcode;
    int counter = 0;
    std::cout << "Loading mappings file " << mapping_filename << std::endl;
    while (std::getline(infile, line)){
        // read name, contig, number of kmers
        std::istringstream(line) >> fields[0] >> fields[1] >> fields[2] ;
        barcode = fields[0].substr(fields[0].find("_") + 1);
        add_barcode_vote(barcode, fields[2], std::stoi(fields[1]));
        counter += 1;

    }
    std::cout << "Loaded " << counter << " mappings from " << barcode_edge_mappings.size() << " barcodes" <<std::endl;
}
