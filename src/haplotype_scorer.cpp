//
// Created by Katie Barr on 11/07/2026.
//

#include <numeric>
#include <unordered_set>
#include "haplotype_scorer.h"



double avg(std::vector<int> const& v){
    if (v.size() > 0) {
        return static_cast<double>(std::accumulate(v.begin(), v.end(), 0LL)) / v.size();
    }
    return 0.0;
}

double stdev(std::vector<int> const& v, double mean){
    if (v.size() > 0) {
        double res = 0;
        for (auto i: v) {
            res += std::pow(i - mean, 2);
        }
        return std::pow(res / v.size(), 0.5);
    }
    return 0.0;
}


HaplotypeScorer::HaplotypeScorer(std::string mapping_file, std::vector<std::vector <std::string> > possible_hs, Graph& g): graph(g){
    mapping_filename=std::move(mapping_file);
    possible_haplotypes=std::move(possible_hs);
    for (auto &hap: possible_haplotypes){
        for (auto &e: hap){
            haplotype_edges.insert(e); // names of all edges appearing in any haplotype
        }
    }

    std::cout << "edges in haps: "<< haplotype_edges.size() <<std::endl;
    std::cout << "mappings " <<  mapping_filename <<std::endl;
    std::cout << std::endl;
}

void print_vector(std::vector<std::string> const& vec){
    for (auto a: vec){
        std::cout << a << " ";
    }
    std::cout << std::endl;
}

void print_int_vector(std::vector<int> const& vec){
    for (auto a: vec){
        std::cout << a << " ";
    }
    std::cout << std::endl;
}

void print_pair_int_vector(std::vector<std::pair<int, int> > const& vec){
    for (auto a: vec){
        std::cout << std::get<0>(a) << " " << std::get<1>(a);
    }
    std::cout << std::endl;
}

struct SupportStats {
    int max_val;
    int min_val;
    double mean;
    double stdev_val;
};

SupportStats compute_stats(const std::vector<int>& v){
    SupportStats s;
    auto minmax = std::minmax_element(v.begin(), v.end());
    s.min_val = *minmax.first;
    s.max_val = *minmax.second;
    s.mean = avg(v);
    s.stdev_val = stdev(v, s.mean);
    return s;
}

void log_stats(const std::string &label, const SupportStats &s){
    std::cout << label << " max : " << s.max_val << " min: " << s.min_val
               << " mean: " << s.mean << " stdev: " << s.stdev_val << std::endl;
}

// writes the top `count` ranked haplotype indices as "<tag>\t<hap_label>\t<rank>\t<edge>\t<edge>..."
void write_ranked_block(std::ofstream &out, const std::string &tag, const std::string &hap_label, size_t count,
                         const std::vector<int> &indices,
                         const std::vector<std::vector<std::string> > &possible_haplotypes){
    for (size_t i = 0; i < count; i++){
        out << tag << "\t" << hap_label << "\t" << i << "\t";
        for (auto h: possible_haplotypes[indices[i]]){
            out << h << "\t";
        }
        out << std::endl;
    }
}

std::vector<int> HaplotypeScorer::winner_for_barcode(const std::unordered_map<int, int>& haplotype_scores){
    int max=0;
    std::vector<int> winners;
    for (auto h: haplotype_scores){
        if (h.second > max){
            max = h.second;
            winners.clear();
            winners.push_back(h.first);
        } else if (h.second == max){
            winners.push_back(h.first);
        }
    }
    return winners;
}

void HaplotypeScorer::print_summary(const std::string& outfile, const std::vector<std::pair<int, int> >& supports, const std::vector<std::pair<int, int> >& overall_supports, const std::vector<int>& haplotype_support_vals, const std::vector<int>& haplotype_overall_support_vals, const std::vector<int>& haplotype_not_support_vals){
    //print results in a command-line parseable way
    // need scores for each, maybe top 3, potentisal hap choices- actual edge seqs for each
    // also include stats
    std::ofstream out(outfile);
    // support score 1, support score 1 hap, overall support score 1, overall support score 1 hap, support score 2, support score 2 hap, overall support score 2, overall support score 2 hap, support score 3, support score 3 hap, overall support score 3, overall support score 3 hap
    for (int i=0; i < 3; i++){
        out << supports[i].first << "\t" << supports[i].second << "\t" << overall_supports[i].first << "\t" <<  overall_supports[i].second;
    }
    out << std::endl;

    SupportStats support_stats = compute_stats(haplotype_support_vals);
    SupportStats not_support_stats = compute_stats(haplotype_not_support_vals);
    SupportStats overall_stats = compute_stats(haplotype_overall_support_vals);

    log_stats("Haplotype support", support_stats);
    log_stats("Haplotype not support", not_support_stats);
    log_stats("Haplotype overall support", overall_stats);

    // max, min, mean, stdev for support, not support, overall support
    out << support_stats.max_val << "\t" << support_stats.min_val << "\t" << support_stats.mean << "\t" << support_stats.stdev_val << "\t"
        << not_support_stats.max_val << "\t" << not_support_stats.min_val << "\t" << not_support_stats.mean << "\t"
        << overall_stats.max_val << "\t" << overall_stats.min_val << std::endl;

    std::vector<int> support_indices, overall_indices;
    for (int i = 0; i < 3; i++){
        support_indices.push_back(supports[i].first);
        overall_indices.push_back(overall_supports[i].first);
    }
    write_ranked_block(out, "Barcode support", "Hap", 3, support_indices, possible_haplotypes);
    write_ranked_block(out, "Overall support", "Hap", 3, overall_indices, possible_haplotypes);

    max_overall_support = support_stats.max_val;
    mean_overall_support = support_stats.mean;
}

// void HaplotypeScorer::print_pair_summary(const std::string& outfile, const std::vector<std::pair<std::pair<int, int>, int > >& pair_supports, const std::vector<std::pair<std::pair<int, int>, int > >& pair_overall_supports, const std::vector<int>& hap_pair_support_values, const std::vector<int>& hap_pair_support_total_score_values, const std::vector<int>& hap_pair_not_support_values){
//     std::ofstream out;
//     out.open( outfile.c_str(),  std::ofstream::out | std::ofstream::app );
//     // exactly as above but for pairs
//     // support score 1, support score 1 hap, overall support score 1, overall support score 1 hap, support score 2, support score 2 hap, overall support score 2, overall support score 2 hap, support score 3, support score 3 hap, overall support score 3, overall support score 3 hap
//     for (int i=0; i < 3; i++){
//         out << std::get<0>(pair_supports[i].first) << "\t" << std::get<1>(pair_supports[i].first) << "\t" << pair_supports[i].second << "\t" <<std::get<0>(pair_overall_supports[i].first) << "\t" << std::get<1>(pair_overall_supports[i].first) << "\t" << pair_overall_supports[i].second;
//     }
//     out << std::endl;

//     SupportStats pair_support_stats = compute_stats(hap_pair_support_values);
//     SupportStats pair_not_support_stats = compute_stats(hap_pair_not_support_values);
//     SupportStats pair_overall_stats = compute_stats(hap_pair_support_total_score_values);

//     std::cout << "pair support size: " << pair_supports.size() << " total: "
//               << hap_pair_support_total_score_values.size() << std::endl;
//     max_overall_pair_support = pair_overall_stats.max_val;
//     mean_overall_pair_support = pair_overall_stats.mean;

//     log_stats("Haplotype pair support", pair_support_stats);
//     log_stats("Haplotype pair not support", pair_not_support_stats);
//     log_stats("Haplotype pair overall support", pair_overall_stats);

//     // max, min, mean, stdev for support, not support, overall support
//     out << pair_support_stats.max_val << "\t" << pair_support_stats.min_val << "\t" << pair_support_stats.mean << "\t" << pair_support_stats.stdev_val << "\t"
//         << pair_not_support_stats.max_val << "\t" << pair_not_support_stats.min_val << "\t" << pair_not_support_stats.mean << "\t" << pair_not_support_stats.stdev_val << "\t"
//         << pair_overall_stats.max_val << "\t" << pair_overall_stats.min_val << std::endl;

//     size_t len = std::min<size_t>(3, pair_supports.size());
//     std::vector<int> p1_support, p2_support, p1_overall, p2_overall;
//     for (size_t i = 0; i < len; i++){
//         p1_support.push_back(std::get<0>(pair_supports[i].first));
//         p2_support.push_back(std::get<1>(pair_supports[i].first));
//         p1_overall.push_back(std::get<0>(pair_overall_supports[i].first));
//         p2_overall.push_back(std::get<1>(pair_overall_supports[i].first));
//     }
//     write_ranked_block(out, "Barcode support", "HapP1", len, p1_support, possible_haplotypes);
//     write_ranked_block(out, "Barcode support", "HapP2", len, p2_support, possible_haplotypes);
//     write_ranked_block(out, "Overall support", "HapP1", len, p1_overall, possible_haplotypes);
//     write_ranked_block(out, "Overall support", "HapP2", len, p2_overall, possible_haplotypes);
// }

int HaplotypeScorer::score_haplotypes(std::string outfile) {
    std::cout << possible_haplotypes.size() << std::endl;
    //initialize score arrays- index is haplotype index
    std::vector<int> haplotype_support(possible_haplotypes.size(), 0);
    std::vector<int> haplotype_not_support(possible_haplotypes.size(), 0);
    // support whether wins or not
    std::vector<int> haplotype_overall_support(possible_haplotypes.size(), 0);
    std::string barcode;
    for (auto &bm: barcode_haplotype_mappings) {
        barcode = bm.first;
        std::vector<int> winners = winner_for_barcode(bm.second); // ideally should be length 1, gives index of possible haplotype with most support
        // winners is a list of indices from possible haplotypes
        for (int i = 0; i < possible_haplotypes.size(); i++) {
            auto score_it = bm.second.find(i);
            int barcode_hap_score = score_it != bm.second.end() ? score_it->second : 0; // score for this barcode for this haplotype
            if (std::find(winners.begin(), winners.end(), i) == winners.end()) {
                haplotype_barcode_disagree[i][barcode] += barcode_hap_score;
                haplotype_not_support[i] += 1;
            }
            else { // if haplotype in winners
                haplotype_support[i] += 1;
                haplotype_barcode_agree[i][barcode] += barcode_hap_score;
            }

            // support whether wins or not
            haplotype_overall_support[i] += barcode_hap_score;
        }

    }

    std::vector<std::pair<int, int> > supports;
    std::vector<std::pair<int, int> > overall_supports;

    for (int i=0; i< possible_haplotypes.size(); i++){
        supports.push_back(std::make_pair(i, haplotype_support[i]));
        overall_supports.push_back(std::make_pair(i, haplotype_overall_support[i]));
    }

    std::sort(supports.begin(), supports.end(), [](auto &left, auto &right) {
        return left.second < right.second;
    });
    std::sort(overall_supports.begin(), overall_supports.end(), [](auto &left, auto &right) {
        return left.second < right.second;
    });
    // when not rushing, get rid of al this repetition
    auto support_max = std::max_element(haplotype_support.begin(), haplotype_support.end());
    auto support_mean = avg(haplotype_support);
    auto overall_support_max = std::max_element(haplotype_overall_support.begin(),
                                                haplotype_overall_support.end());
    auto overall_support_mean = avg(haplotype_overall_support);

    auto support_stdev = stdev(haplotype_support, support_mean);

    auto overall_stdev = stdev(haplotype_overall_support, overall_support_mean);

    print_summary(outfile, supports, overall_supports, haplotype_support, haplotype_overall_support, haplotype_not_support);

    // get winners
    std::vector<int> support_winner;
    std::vector<int> overall_support_winner;
    for (int h = 0; h < possible_haplotypes.size(); h++) {
        if (haplotype_support[h] == *support_max) {
            support_winner.push_back(h);
        }
        if (haplotype_overall_support[h] == *overall_support_max) {
            overall_support_winner.push_back(h);
        }
    }

    std::cout << "Support winner: ";
    print_int_vector(support_winner);
    std::cout << "overall SUpport winner: ";
    print_int_vector(overall_support_winner);
    print_vector(possible_haplotypes[support_winner[0]]);

    if (*support_max > *overall_support_max) { // if the support winner is also the overall support winner, we have a clear winner
        winners = std::make_pair(possible_haplotypes[support_winner[0]], possible_haplotypes[overall_support_winner[0]]);
        winning_haplotype = std::make_pair(support_winner[0], overall_support_winner[0]);
        return 0;
    } else  if (*overall_support_max > (overall_support_mean + 2 *overall_stdev)) { // if it doesn't make it... pick best we can do, so overall support
        winning_haplotype = std::make_pair(support_winner[0], overall_support_winner[0]);
        return 1;
    } else {
        return 2;
    }
return 2;
}


void HaplotypeScorer::decide_barcode_haplotype_support(){

    std::vector<std::unordered_set<std::string> > haplotype_edge_sets(possible_haplotypes.size());
    for (size_t i = 0; i < possible_haplotypes.size(); i++){
        haplotype_edge_sets[i] = {possible_haplotypes[i].begin(), possible_haplotypes[i].end()};
    }
    std::cout << "Calculating barcode haplotype support for " << barcode_edge_mappings.size() << " mappings"<< std::endl;
    for (auto &mapping:barcode_edge_mappings){
        // if barcode maps to more than 1 edge in bubbles and maximum support is greater than 1
        if (mapping.second.size() > 1){ // if barcode maps to more than one edge
            int max_score = 0;
            for (auto &e: mapping.second){
                max_score = std::max(max_score, e.second);
            }
            if (max_score > 1) { // is score > 1
                for (int i = 0; i < possible_haplotypes.size(); i++) {
                    auto& hset = haplotype_edge_sets[i];
                    // count edges (and sum their support) that this barcode maps to and that belong to this haplotype
                    size_t match_count = 0;
                    int support = 0;
                    for (auto &e: mapping.second) {
                        if (hset.count(e.first)) {
                            match_count += 1;
                            support += e.second;
                        }
                    }
                    // somewhat arbitrary rule to decide if the barcode supports a haplotype enough
                    if (match_count >= (mapping.second.size() / 2) && match_count > 1) {
                        barcode_haplotype_mappings[mapping.first][i] = support;
                    }
                }
            }

        } else {
            unused_barcodes.push_back(mapping.first);
        }
    }
    std::cout << "Calculated haplotype support for each barcode, " << barcode_haplotype_mappings.size() <<  std::endl;

}

void HaplotypeScorer::write_output_partial_success(std::string output_file){
    std::string o = "partial_" + output_file;// output is same for partial, just need to know somehow that were less confident
    write_output_success(o);
}
void HaplotypeScorer::write_output_success(std::string output_file){
    std::ofstream out(output_file + ".txt");
    std::vector<std::string> winner1 = std::get<0>(winners);
    std::vector<std::string> winner2 = std::get<1>(winners);
    out << "Haplotype 1: " << std::endl;
    for (auto h: winner1){
        out << h << " ";
    }
    out << std::endl;
    out << "Haplotype 2: " << std::endl;
    for (auto h: winner2){
        out << h << " ";
    }
    graph.output_contigs_joined_to_contig_list(winner1, haplotype_barcode_agree[std::get<0>(winning_haplotype)], output_file + ".hap1");
    graph.output_contigs_joined_to_contig_list(winner2, haplotype_barcode_agree[std::get<1>(winning_haplotype)], output_file + ".hap2");
    out << std::endl;
    out << "Overall support for pair: " << max_overall_pair_support << " mean:" << mean_overall_pair_support <<std::endl;
    out << "Highest overall individual hap support: " << max_overall_support << " mean: " << mean_overall_support <<std::endl;
    out << "Barcodes supporting winner, hap1:" << std::endl;
    std::unordered_set<std::string> barcodes_seen;
    // need barcodes supporting this pair- to outputm for each barcode, total kmers, kmers agreeing, kmers disagreeing, kmers to hom parts, other
    for (auto &b:haplotype_barcode_agree[std::get<0>(winning_haplotype)]){
        barcodes_seen.insert(b.first);
        int total_agreeing_kmers = b.second;
        int total_hom_kmers = barcode_hom_mappings[b.first];
        int total_kmers = kmers_per_barcode[b.first];
        int total_disagreeing_kmers = haplotype_barcode_disagree[std::get<0>(winning_haplotype)][b.first];
        // kmers mapping elsewhere is just total minus all others
        int other = total_kmers - total_agreeing_kmers - total_hom_kmers - total_disagreeing_kmers;
        out << b.first << ": " << total_agreeing_kmers << ", " << total_disagreeing_kmers << ", " << total_hom_kmers <<", " << other << ", " << total_kmers << std::endl;
    }
    out << "Barcodes supporting winner, hap2:" << std::endl;
    // need barcodes supporting this pair- to outputm for each barcode, total kmers, kmers agreeing, kmers disagreeing, kmers to hom parts, other
    for (auto &b:haplotype_barcode_agree[std::get<1>(winning_haplotype)]){
        barcodes_seen.insert(b.first);
        int total_agreeing_kmers = b.second;
        int total_hom_kmers = barcode_hom_mappings[b.first];
        int total_kmers = kmers_per_barcode[b.first];
        int total_disagreeing_kmers = haplotype_barcode_disagree[std::get<1>(winning_haplotype)][b.first];
        // kmers mapping elsewhere is just total minus all others
        int other = total_kmers - total_agreeing_kmers - total_hom_kmers - total_disagreeing_kmers;
        out << b.first << ": " << total_agreeing_kmers << ", " << total_disagreeing_kmers << ", " << total_hom_kmers <<", " << other << ", " << total_kmers << std::endl;
    }
    // then need other barcodes which mapped usefully to this region but didn't support
    out << "Barcodes mapping to this region that do not support winner:" << std::endl;
    for (auto &b: barcode_haplotype_mappings){
        if (barcodes_seen.find(b.first) == barcodes_seen.end()){
            int total_agreeing_kmers = barcode_haplotype_mappings[b.first][std::get<0>(winning_haplotype)] + barcode_haplotype_mappings[b.first][std::get<1>(winning_haplotype)];
            int total_hom_kmers = barcode_hom_mappings[b.first];
            int total_kmers = kmers_per_barcode[b.first];
            int total_disagreeing_kmers = haplotype_barcode_disagree[std::get<1>(winning_haplotype)][b.first];
            int other = total_kmers - total_agreeing_kmers - total_hom_kmers - total_disagreeing_kmers;
            out << b.first << ": " << total_agreeing_kmers << ", " << total_disagreeing_kmers << ", " << total_hom_kmers <<", " << other << ", " << total_kmers << std::endl;
        }
    }
}


void HaplotypeScorer::add_barcode_vote(const std::string& barcode, const std::string& edge, int kmers){
    barcodes.insert(barcode);
    if (haplotype_edges.find(edge) != haplotype_edges.end()){
        // we only care about mappings to edges in bubbles, which will all be in haplotype_edges
        barcode_edge_mappings[barcode][edge] += kmers;
    } else if (graph.edges.find(edge) != graph.edges.end()){
        barcode_hom_mappings[barcode] += kmers;
    }
    kmers_per_barcode[barcode] += kmers;

}

void HaplotypeScorer::load_mappings_from_dict(std::map<std::string, std::map<std::string, int> > & mappings) {
    for (auto &edge:graph.edges){
        auto it = mappings.find(edge);
        if (it == mappings.end()) continue;
        for (auto &barcode: it->second){
            add_barcode_vote(barcode.first, edge, barcode.second);
        }
    }
}
