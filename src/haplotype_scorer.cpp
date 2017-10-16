//
// Created by Katie Barr (EI) on 15/10/2017.
//

#include <numeric>
#include "haplotype_scorer.h"
double avg(std::vector<int> v){
    return std::accumulate(v.begin(), v.end(), 0LL) / v.size();
}

double stdev(std::vector<int> v, double mean){
    double res=0;
    for (auto i: v){
        res += std::pow(i-mean,2);
    }
    return std::pow(res/v.size(), 0.5);
}


HaplotypeScorer::HaplotypeScorer(std::string mapping_file, std::vector<std::vector <std::string> > possible_hs, Graph g){
    mapping_filename=mapping_file;
    possible_haplotypes=possible_hs;
    graph = g;
    std::set<std::string> edges;
    for (auto hap: possible_haplotypes){
        for (auto e: hap){
            edges.insert(e);
        }
    }
    std::cout << "edges in haps: "<< edges.size() <<std::endl;
    for (auto e:edges){
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
        std::cout << a << " ";
    }
    std::cout << std::endl;
}

void print_int_vector(std::vector<int> vec){
    for (auto a: vec){
        std::cout << a << " ";
    }
    std::cout << std::endl;
}

void print_pair_int_vector(std::vector<std::pair<int, int> > vec){
    for (auto a: vec){
        std::cout << std::get<0>(a) << " " << std::get<1>(a);
    }
    std::cout << std::endl;
}

std::vector<int>  HaplotypeScorer::winner_for_barcode(std::string barcode){
    int max=0;
    std::vector<int> winners;
    for (auto h:barcode_haplotype_mappings[barcode]){
        if (h.second > max){
            max = h.second;
        }
    }
    //TODO: DECIDE CRITERIA FOR MINIMUM SUPPORT
    if (max <10){
        return winners;
    }
    for (auto h:barcode_haplotype_mappings[barcode]){
        if (h.second == max){
            winners.push_back(h.first);
        }
    }
    return winners;
}

std::vector<std::string> HaplotypeScorer::score_haplotypes(){
    //initialize score arrays- index is haplotype index
    int haplotype_support[possible_haplotypes.size()] = {0};
    int haplotype_not_support[possible_haplotypes.size()] = {0};
    int haplotype_overall_support[possible_haplotypes.size()] = {0};
    std::map<std::pair<int, int>, int> hap_pair_not_support;
    std::map<std::pair<int, int>, int>  hap_pair_support;
    std::map<std::pair<int, int>, int> hap_pair_support_total_score;
    std::string barcode;
    std::vector<int> winners;
    int pair;
    for (auto bm: barcode_haplotype_mappings){
        barcode = bm.first;
        winners = winner_for_barcode(barcode);
        if (winners.size() > 0) {
            for (int hap = 0; hap < possible_haplotypes.size() / 2; hap++) {
                pair = possible_haplotypes.size() / 2;
                if (bm.second.find(hap) != bm.second.end() or bm.second.find(pair) != bm.second.end() ) {
                    if (bm.second.find(hap) != bm.second.end() ){
                        haplotype_overall_support[hap] += bm.second[hap];
                        if (std::find(winners.begin(), winners.end(), hap) != winners.end()){
                            haplotype_support[hap] += 1;
                            hap_pair_support[std::make_pair(hap, pair)] += 1;
                            hap_pair_support_total_score[std::make_pair(hap, pair)] += bm.second[hap];
                        } else {
                            haplotype_not_support[hap] += 1;
                        }

                    } else if (bm.second.find(pair) != bm.second.end() ){
                        haplotype_overall_support[pair] += bm.second[pair];
                        if (std::find(winners.begin(), winners.end(), pair) != winners.end()){
                            haplotype_support[pair] += 1;
                            hap_pair_support[std::make_pair(hap, pair)] += 1;
                            hap_pair_support_total_score[std::make_pair(hap, pair)] += bm.second[pair];
                        } else {
                            haplotype_not_support[hap] += 1;
                        }
                    }
                } else if (bm.second.find(hap) == bm.second.end() and bm.second.find(pair) == bm.second.end() ){
                    hap_pair_not_support[std::make_pair(hap, pair)] += 1;
                }
            }
        }
    }
    std::vector<int > haplotype_support_vals;
    std::vector<int > haplotype_not_support_vals;
    std::vector<int > haplotype_overall_support_vals;
    for (int i=0; i < possible_haplotypes.size(); i++){
        haplotype_support_vals.push_back(haplotype_support[i]);
        haplotype_not_support_vals.push_back(haplotype_not_support[i]);
        haplotype_overall_support_vals.push_back(haplotype_overall_support[i]);

    }
    auto support_max = std::max_element(haplotype_support_vals.begin(), haplotype_support_vals.end());
    auto support_mean = avg(haplotype_support_vals);
    auto overall_support_max = std::max_element(haplotype_overall_support_vals.begin(), haplotype_overall_support_vals.end());
    auto overall_support_mean = avg(haplotype_support_vals);

    auto support_stdev = stdev(haplotype_support_vals, support_mean);

    auto overall_stdev = stdev(haplotype_overall_support_vals, overall_support_mean);
    std::cout << "Haplotype support max : " << *support_max << " min: " << *std::min_element(haplotype_support_vals.begin(), haplotype_support_vals.end()) << " mean: "<< support_mean << " stdev: " << support_stdev << std::endl;
    std::cout << "Haplotype not support max : " << *std::max_element(haplotype_not_support_vals.begin(), haplotype_not_support_vals.end()) << " min: " << *std::min_element(haplotype_not_support_vals.begin(), haplotype_not_support_vals.end()) << " mean: "<< avg(haplotype_not_support_vals)<< std::endl;
    std::cout << "Haplotype overall support max: " << *overall_support_max << " min: " << *std::min_element(haplotype_overall_support_vals.begin(), haplotype_overall_support_vals.end()) << " mean: "<< overall_support_mean << " stdev: " << overall_stdev << std::endl;

    std::vector< int> hap_pair_support_values;
    for (auto h : hap_pair_support){
        hap_pair_support_values.push_back(h.second);
    }
    std::vector< int> hap_pair_not_support_values;
    for (auto h : hap_pair_not_support){
        hap_pair_not_support_values.push_back(h.second);
    }
    std::vector< int> hap_pair_support_total_score_values;
    for (auto h : hap_pair_support_total_score){
        hap_pair_support_total_score_values.push_back(h.second);
    }
    auto pair_support_max = std::max_element(hap_pair_support_values.begin(), hap_pair_support_values.end());
    auto not_pair_support_max = std::max_element(hap_pair_not_support_values.begin(), hap_pair_not_support_values.end());
    auto overall_pair_support_max = std::max_element(hap_pair_support_total_score_values.begin(), hap_pair_support_total_score_values.end());
    auto pair_support_mean = avg(hap_pair_support_values);
    auto pair_support_stdev = stdev(hap_pair_support_values, pair_support_mean);
    auto pair_not_support_mean = avg(hap_pair_not_support_values);
    auto pair_not_support_stdev = stdev(hap_pair_not_support_values, pair_not_support_mean);
    auto pair_overall_support_mean = avg(hap_pair_support_values);
    auto pair_overall_support_stdev = stdev(hap_pair_support_total_score_values, pair_overall_support_mean);
    std::cout << "Haplotype pair support max : " << *pair_support_max << " mean: " << pair_support_mean << " stdev: " << pair_support_stdev<< std::endl;
    std::cout << "Haplotype pair not support max : " << *not_pair_support_max << " mean: " << pair_not_support_mean << " stdev: " << pair_not_support_stdev<< std::endl;
    std::cout << "Haplotype pair overall support max : " << *overall_pair_support_max << " mean: " << pair_overall_support_mean << " stdev: " << pair_overall_support_stdev<< std::endl;

    // get winners
    std::vector<int> support_winner;
    std::vector<int>  overall_support_winner;
    for (int h=0; h < possible_haplotypes.size(); h++){
        if (haplotype_support[h] == *support_max){
            support_winner.push_back(h);
        }
        if (haplotype_overall_support[h] == *overall_support_max){
            overall_support_winner.push_back(h);
        }
    }
    std::vector<std::pair<int, int> > pair_support_winner;
    std::vector<std::pair<int, int> > pair_overall_support_winner;
    for (auto h: hap_pair_support){
        if (h.second == *pair_support_max){
            pair_support_winner.push_back(h.first);
        }
        if (h.second == *overall_pair_support_max){
            pair_overall_support_winner.push_back(h.first);
        }
    }
    std::cout << "Support winner: ";
    print_int_vector(support_winner);
    std::cout << "overall SUpport winner: ";
    print_int_vector(overall_support_winner);
    std::cout << "pair SUpport winner: ";
    print_pair_int_vector(pair_support_winner);
    std::cout << "pair overall SUpport winner: ";
    print_pair_int_vector(pair_overall_support_winner);
    print_vector(possible_haplotypes[support_winner[0]]);


}


void HaplotypeScorer::decide_barcode_haplotype_support(){
    std::vector<std::string> edges;
    std::vector<std::string> edges_in_haplotype;
    std::vector<std::string> h;
    int support;
    int haplotypes_supported = 0;
    std::cout << "Calculating barcode haplotype support" << std::endl;
    for (auto mapping:barcode_edge_mappings){
        //std::cout << "Checking barcode " << mapping.first <<std::endl;
        // if barcode maps to more than 1 edge in bubbles
        if (mapping.second.size() > 1){
            for (auto e: mapping.second){
                edges.push_back(e.first);
            }
            for (int i=0; i < possible_haplotypes.size(); i++){
                h = possible_haplotypes[i];
                // find all edges in each haplotype that this barcode maps to
                for (auto e1: h){
                    if (std::find(edges.begin(), edges.end(), e1) != edges.end()){
                        edges_in_haplotype.push_back(e1);
                    }
                }

                // somewhat arbitrary rule to decide if the barcode supports a haplotype enough
                if (edges_in_haplotype.size() > edges.size()/2 && edges_in_haplotype.size() >1){
                    for (auto a: edges_in_haplotype){
                        support += mapping.second[a];
                    }
                    barcode_haplotype_mappings[mapping.first][i] = support;
                    support = 0;
                    haplotypes_supported +=1;
                } else {
                    unused_barcodes.push_back(mapping.first);
                }
                edges_in_haplotype.clear();
            }

            edges.clear();

        } else {
            unused_barcodes.push_back(mapping.first);
        }
        //std::cout << "barcode " << mapping.first << " supports " << haplotypes_supported << std::endl;

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
