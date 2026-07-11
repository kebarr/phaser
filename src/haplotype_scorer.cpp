//
// Created by Katie Barr on 11/07/2026.
//

#include <numeric>
#include <unordered_set>
#include "haplotype_scorer.h"

template <typename T, typename T2=T>
struct accumulator
{
    T2 sum; // we could plug in a more accurate type for the sum
    T S;
    T M;
    size_t N;

    // default constructor initializes all values
    accumulator() : sum(0), S(0), M(0), N(0) { }

    // add another number
    T2 operator()(const T& x) {
        ++N;
        sum += x;
        T Mprev = M;
        M += (x - Mprev) / N;
        S += (x - Mprev) * (x - M);
        return sum;
    }

    T mean() const { return sum / N; }

    T variance() const { return S / (N - 1); }

    // operator<< to print the statistics to screen:
    // denoted friend just to be able to write this inside
    // the class definition and thus not to need to write
    // the template specification of accumulator...
    friend std::ostream& operator<<(std::ostream& out,
                                    const accumulator& a)
    {
        if (a.N > 0)
            out << "N\t\t\t= " << a.N << std::endl
                << "sum\t\t\t= " << a.sum << std::endl
                << "mean\t\t= " << std::fixed << std::setprecision(2) << a.mean() << std::endl;
        if (a.N > 1)
            out << "sd\t\t\t= " << std::fixed << std::setprecision(2) << std::sqrt(a.variance()) << std::endl;
        else
            out << "sd\t\t\t= " << std::fixed << std::setprecision(2) << 0 << std::endl;
        return out;
    }

};


double avg(std::vector<int> v){
    if (v.size() > 0) {
        return std::accumulate(v.begin(), v.end(), 0LL) / v.size();
    }
    return 0.0;
}

double stdev(std::vector<int> v, double mean){
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
    mapping_filename=mapping_file;
    possible_haplotypes=possible_hs;
    std::set<std::string> edges;
    for (auto hap: possible_haplotypes){
        for (auto e: hap){
            edges.insert(e); // names of all edges
        }
    }

    std::cout << "edges in haps: "<< edges.size() <<std::endl;
    for (auto e:edges){
        for (int i=0; i < possible_haplotypes.size(); i++){
                if (std::find(possible_haplotypes[i].begin(), possible_haplotypes[i].end(), e) != possible_haplotypes[i].end()){
                    edge_haplotype_dict[e].push_back(i); /// if e is in this possible haplotype then add its index
                }
        }
    }
    std::cout << "mappings " <<  mapping_filename <<std::endl;
    std::cout << "edges in haplotypes: " << edge_haplotype_dict.size() <<std::endl;
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

struct SupportStats {
    int max_val;
    int min_val;
    double mean;
    double stdev_val;
};

SupportStats compute_stats(const std::vector<int> &v){
    SupportStats s;
    s.max_val = *std::max_element(v.begin(), v.end());
    s.min_val = *std::min_element(v.begin(), v.end());
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

std::vector<int>  HaplotypeScorer::winner_for_barcode(std::string barcode){
    int max=0;
    std::vector<int> winners;
    for (auto h:barcode_haplotype_mappings[barcode]){
        if (h.second > max){
            max = h.second;
        }
    }
    //TODO: DECIDE CRITERIA FOR MINIMUM SUPPORT
    /*if (max <10){
        return winners;
    }*/
    for (auto h:barcode_haplotype_mappings[barcode]){
        if (h.second == max){
            winners.push_back(h.first);
        }
    }
    return winners;
}

void HaplotypeScorer::print_summary(std::string outfile, std::vector<std::pair<int, int> > supports, std::vector<std::pair<int, int> > overall_supports, std::vector<int> haplotype_support_vals, std::vector<int> haplotype_overall_support_vals, std::vector<int> haplotype_not_support_vals){
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

void HaplotypeScorer::print_pair_summary(std::string outfile, std::vector<std::pair<std::pair<int, int>, int > > pair_supports, std::vector<std::pair<std::pair<int, int>, int > > pair_overall_supports, std::vector<int> hap_pair_support_values, std::vector<int> hap_pair_support_total_score_values, std::vector<int> hap_pair_not_support_values){
    std::ofstream out;
    out.open( outfile.c_str(),  std::ofstream::out | std::ofstream::app );
    // exactly as above but for pairs
    // support score 1, support score 1 hap, overall support score 1, overall support score 1 hap, support score 2, support score 2 hap, overall support score 2, overall support score 2 hap, support score 3, support score 3 hap, overall support score 3, overall support score 3 hap
    for (int i=0; i < 3; i++){
        out << std::get<0>(pair_supports[i].first) << "\t" << std::get<1>(pair_supports[i].first) << "\t" << pair_supports[i].second << "\t" <<std::get<0>(pair_overall_supports[i].first) << "\t" << std::get<1>(pair_overall_supports[i].first) << "\t" << pair_overall_supports[i].second;
    }
    out << std::endl;

    SupportStats pair_support_stats = compute_stats(hap_pair_support_values);
    SupportStats pair_not_support_stats = compute_stats(hap_pair_not_support_values);
    SupportStats pair_overall_stats = compute_stats(hap_pair_support_total_score_values);

    std::cout << "pair support size: " << pair_supports.size() << " total: "
              << hap_pair_support_total_score_values.size() << std::endl;
    max_overall_pair_support = pair_overall_stats.max_val;
    mean_overall_pair_support = pair_overall_stats.mean;

    log_stats("Haplotype pair support", pair_support_stats);
    log_stats("Haplotype pair not support", pair_not_support_stats);
    log_stats("Haplotype pair overall support", pair_overall_stats);

    // max, min, mean, stdev for support, not support, overall support
    out << pair_support_stats.max_val << "\t" << pair_support_stats.min_val << "\t" << pair_support_stats.mean << "\t" << pair_support_stats.stdev_val << "\t"
        << pair_not_support_stats.max_val << "\t" << pair_not_support_stats.min_val << "\t" << pair_not_support_stats.mean << "\t" << pair_not_support_stats.stdev_val << "\t"
        << pair_overall_stats.max_val << "\t" << pair_overall_stats.min_val << std::endl;

    size_t len = std::min<size_t>(3, pair_supports.size());
    std::vector<int> p1_support, p2_support, p1_overall, p2_overall;
    for (size_t i = 0; i < len; i++){
        p1_support.push_back(std::get<0>(pair_supports[i].first));
        p2_support.push_back(std::get<1>(pair_supports[i].first));
        p1_overall.push_back(std::get<0>(pair_overall_supports[i].first));
        p2_overall.push_back(std::get<1>(pair_overall_supports[i].first));
    }
    write_ranked_block(out, "Barcode support", "HapP1", len, p1_support, possible_haplotypes);
    write_ranked_block(out, "Barcode support", "HapP2", len, p2_support, possible_haplotypes);
    write_ranked_block(out, "Overall support", "HapP1", len, p1_overall, possible_haplotypes);
    write_ranked_block(out, "Overall support", "HapP2", len, p2_overall, possible_haplotypes);
}

int HaplotypeScorer::score_haplotypes(std::string outfile) {
    std::cout << possible_haplotypes.size() << std::endl;
    //initialize score arrays- index is haplotype index
    int haplotype_support[possible_haplotypes.size()] = {0};
    int haplotype_not_support[possible_haplotypes.size()] = {0};
    int haplotype_overall_support[possible_haplotypes.size()] = {0};
    std::map<std::pair<int, int>, int> hap_pair_not_support;
    std::map<std::pair<int, int>, int> hap_pair_support;
    std::map<std::pair<int, int>, int> hap_pair_support_total_score;
    std::string barcode;
    for (auto &bm: barcode_haplotype_mappings) {
        barcode = bm.first;
        std::vector<int> winners = winner_for_barcode(barcode); // ideally should be length 1
        for (auto winner:winners){
            int pair = possible_haplotypes.size() - 1 - winner;
            haplotype_support[winner] += 1;
            hap_pair_support[std::make_pair(winner, pair)] += 1;
            haplotype_barcode_agree[winner][barcode] += bm.second[winner];
            haplotype_barcode_disagree[winner][barcode] += bm.second[pair];
        }
        // for haplotype in range(len(self.list_of_possible_haplotypes)/2):
        for (int hap = 0; hap < possible_haplotypes.size() / 2; hap++) {
            // pair = len(self.list_of_possible_haplotypes) -1 -haplotype
            int pair = possible_haplotypes.size() - 1 - hap;
            if (bm.second.find(hap) != bm.second.end()) {
                haplotype_overall_support[hap] += bm.second[hap];
                hap_pair_support_total_score[std::make_pair(hap, pair)] += bm.second[hap];
            }

            if (bm.second.find(pair) != bm.second.end()) {
                haplotype_overall_support[pair] += bm.second[pair];
                hap_pair_support_total_score[std::make_pair(hap, pair)] += bm.second[pair];
            }
            if (bm.second.find(hap) == bm.second.end()) {
                haplotype_not_support[hap] += 1;
            }
            if (bm.second.find(pair) == bm.second.end()) {
                haplotype_not_support[pair] += 1;
            }
            if (bm.second.find(hap) == bm.second.end() and bm.second.find(pair) == bm.second.end()) {
                hap_pair_not_support[std::make_pair(hap, pair)] += 1;

            }
        }
    }
    std::vector<int> haplotype_support_vals;
    std::vector<int> haplotype_not_support_vals;
    std::vector<int> haplotype_overall_support_vals;
    for (int i = 0; i < possible_haplotypes.size(); i++) {
        haplotype_support_vals.push_back(haplotype_support[i]);
        haplotype_not_support_vals.push_back(haplotype_not_support[i]);
        haplotype_overall_support_vals.push_back(haplotype_overall_support[i]);

    }
    std::vector<std::pair<int, int> > supports;
    std::vector<std::pair<int, int> > overall_supports;
    std::vector<std::pair<std::pair<int, int>, int > > pair_supports;
    std::vector<std::pair<std::pair<int, int>, int > > pair_overall_supports;
    for (int i=0; i< possible_haplotypes.size(); i++){
        supports.push_back(std::make_pair(i, haplotype_support[i]));
        overall_supports.push_back(std::make_pair(i, haplotype_support[i]));
    }
    std::sort(supports.begin(), supports.end(), [](auto &left, auto &right) {
        return left.second < right.second;
    });
    std::sort(overall_supports.begin(), overall_supports.end(), [](auto &left, auto &right) {
        return left.second < right.second;
    });
    // when not rushing, get rid of al this repetition
    auto support_max = std::max_element(haplotype_support_vals.begin(), haplotype_support_vals.end());
    auto support_mean = avg(haplotype_support_vals);
    auto overall_support_max = std::max_element(haplotype_overall_support_vals.begin(),
                                                haplotype_overall_support_vals.end());
    auto overall_support_mean = avg(haplotype_support_vals);

    auto support_stdev = stdev(haplotype_support_vals, support_mean);

    auto overall_stdev = stdev(haplotype_overall_support_vals, overall_support_mean);

    print_summary(outfile, supports, overall_supports, haplotype_support_vals, haplotype_overall_support_vals, haplotype_not_support_vals);
    if (hap_pair_support.size() > 0 && hap_pair_support_total_score.size() > 0) {
        std::vector<int> hap_pair_not_support_values;
        std::vector<int> hap_pair_support_values;
        for (auto h : hap_pair_support) {
            hap_pair_support_values.push_back(h.second);
            pair_supports.push_back(std::make_pair(h.first, h.second));
        }
        for (auto h : hap_pair_not_support) {
            hap_pair_not_support_values.push_back(h.second);
        }
        std::vector<int> hap_pair_support_total_score_values;
        for (auto h : hap_pair_support_total_score) {
            hap_pair_support_total_score_values.push_back(h.second);
            pair_overall_supports.push_back(std::make_pair(h.first, h.second));
        }
        std::sort(pair_supports.begin(), pair_supports.end(), [](auto &left, auto &right) {
            return left.second < right.second;
        });
        std::sort(pair_overall_supports.begin(), pair_overall_supports.end(), [](auto &left, auto &right) {
            return left.second < right.second;
        });

        print_pair_summary(outfile, pair_supports, pair_overall_supports, hap_pair_support_values, hap_pair_support_total_score_values, hap_pair_not_support_values);

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
        std::vector<std::pair<int, int> > pair_support_winner;
        std::vector<std::pair<int, int> > pair_overall_support_winner;
        auto overall_pair_support_max = std::max_element(hap_pair_support_total_score_values.begin(),hap_pair_support_total_score_values.end());

        auto pair_support_max = std::max_element(hap_pair_support_values.begin(), hap_pair_support_values.end());

        for (auto h: hap_pair_support) {
            if (h.second == *pair_support_max) {
                pair_support_winner.push_back(h.first);
            }
        }
        for (auto h: hap_pair_support_total_score) {
            if (h.second == *overall_pair_support_max) {
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

        auto pair_overall_support_mean = avg(hap_pair_support_total_score_values);
        auto pair_overall_support_stdev = stdev(hap_pair_support_total_score_values, pair_overall_support_mean);
        // if they agree on all scores, call it
        if ((std::get<0>(pair_overall_support_winner[0]) == overall_support_winner[0] ||
             std::get<1>(pair_overall_support_winner[0]) == overall_support_winner[0]) &&
            (std::get<0>(pair_support_winner[0]) == support_winner[0] ||
             std::get<1>(pair_support_winner[0]) == support_winner[0])) {
            winners = std::make_pair(possible_haplotypes[std::get<1>(pair_overall_support_winner[0])],
                                     possible_haplotypes[std::get<0>(pair_overall_support_winner[0])]);
            winning_pair = pair_overall_support_winner[0];
            return 0;
        } else if (*overall_pair_support_max > (pair_overall_support_mean + 2 *
                                                                            pair_overall_support_stdev)) { // if it doesn't make it... pick best we can do, so pair overall support
            winners = std::make_pair(possible_haplotypes[std::get<1>(pair_overall_support_winner[0])],
                                     possible_haplotypes[std::get<0>(pair_overall_support_winner[0])]);
            winning_pair = pair_overall_support_winner[0];
            return 1;
        } else {
            return 2;
        }
    }
    return 2;
}


void HaplotypeScorer::decide_barcode_haplotype_support(){

    int support;
    int haplotypes_supported = 0;
    std::vector<std::unordered_set<std::string> > haplotype_edge_sets(possible_haplotypes.size());
    for (size_t i = 0; i < possible_haplotypes.size(); i++){
        haplotype_edge_sets[i] = {possible_haplotypes[i].begin(), possible_haplotypes[i].end()};
    }
    std::cout << "Calculating barcode haplotype support for " << barcode_edge_mappings.size() << " mappings"<< std::endl;
    for (auto &mapping:barcode_edge_mappings){
        //std::cout << "Checking barcode " << mapping.first <<std::endl;
        // if barcode maps to more than 1 edge in bubbles and maximum support is greater than 1
        //auto edge_support_max = std::max_element(std::begin(mapping.second), std::end(mapping.second), [] ( std::map<std::string, int> &p1,  std::map<std::string, int> &p2) {return p1.second < p2.second});
        if (mapping.second.size() > 1){ // if barcode maps to more than one edge
            std::vector<std::string> edges;
            std::vector<int> scores;
            for (auto e: mapping.second){ // add each edge that barcode traverses and its score
                edges.push_back(e.first); 
                scores.push_back(e.second);
            }
            if (*std::max_element(scores.begin(), scores.end())> 1) { 
                for (int i = 0; i < possible_haplotypes.size(); i++) {
                    auto& hset = haplotype_edge_sets[i];
                    // find all edges in each haplotype that this barcode maps to
                    std::vector<std::string> edges_in_haplotype;
                    std::copy_if(edges.begin(), edges.end(), std::back_inserter(edges_in_haplotype),
                                 [&hset](const std::string& e1){ return hset.count(e1) > 0; });
                    // somewhat arbitrary rule to decide if the barcode supports a haplotype enough
                    if (edges_in_haplotype.size() >= (edges.size() / 2) && edges_in_haplotype.size() > 1) {
                        support = 0;
                        for (auto a: edges_in_haplotype) {
                            support += mapping.second[a];
                        }
                        barcode_haplotype_mappings[mapping.first][i] = support;
                        support = 0;
                        haplotypes_supported += 1;
                    } else {
                        unused_barcodes.push_back(mapping.first);
                    }
                }
            }

        } else {
            unused_barcodes.push_back(mapping.first);
        }
        haplotypes_supported = 0;
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
    graph.output_contigs_joined_to_contig_list(winner1, haplotype_barcode_agree[std::get<0>(winning_pair)], output_file + ".hap1");
    graph.output_contigs_joined_to_contig_list(winner2, haplotype_barcode_agree[std::get<1>(winning_pair)], output_file + ".hap2");
    out << std::endl;
    out << "Overall support for pair: " << max_overall_pair_support << " mean:" << mean_overall_pair_support <<std::endl;
    out << "Highest overall individual hap support: " << max_overall_support << " mean: " << mean_overall_support <<std::endl;
    out << "Barcodes supporting winner, hap1:" << std::endl;
    std::vector<std::string> barcodes_seen;
    // need barcodes supporting this pair- to outputm for each barcode, total kmers, kmers agreeing, kmers disagreeing, kmers to hom parts, other
    for (auto b:haplotype_barcode_agree[std::get<0>(winning_pair)]){
        barcodes_seen.push_back(b.first);
        int total_agreeing_kmers = b.second;
        int total_hom_kmers = barcode_hom_mappings[b.first];
        int total_kmers = kmers_per_barcode[b.first];
        int total_disagreeing_kmers = haplotype_barcode_disagree[std::get<0>(winning_pair)][b.first];
        // kmers mapping elsewhere is just total minus all others
        int other = total_kmers - total_agreeing_kmers - total_hom_kmers - total_disagreeing_kmers;
        out << b.first << ": " << total_agreeing_kmers << ", " << total_disagreeing_kmers << ", " << total_hom_kmers <<", " << other << ", " << total_kmers << std::endl;
    }
    out << "Barcodes supporting winner, hap2:" << std::endl;
    // need barcodes supporting this pair- to outputm for each barcode, total kmers, kmers agreeing, kmers disagreeing, kmers to hom parts, other
    for (auto b:haplotype_barcode_agree[std::get<1>(winning_pair)]){
        barcodes_seen.push_back(b.first);
        int total_agreeing_kmers = b.second;
        int total_hom_kmers = barcode_hom_mappings[b.first];
        int total_kmers = kmers_per_barcode[b.first];
        int total_disagreeing_kmers = haplotype_barcode_disagree[std::get<1>(winning_pair)][b.first];
        // kmers mapping elsewhere is just total minus all others
        int other = total_kmers - total_agreeing_kmers - total_hom_kmers - total_disagreeing_kmers;
        out << b.first << ": " << total_agreeing_kmers << ", " << total_disagreeing_kmers << ", " << total_hom_kmers <<", " << other << ", " << total_kmers << std::endl;
    }
    // then need other barcodes which mapped usefully to this region but didn't support
    out << "Barcodes mapping to this region that do not support winner:" << std::endl;
    for (auto b: barcode_haplotype_mappings){
        if (std::find(barcodes_seen.begin(), barcodes_seen.end(), b.first) == barcodes_seen.end()){
            int total_agreeing_kmers = barcode_haplotype_mappings[b.first][std::get<0>(winning_pair)] + barcode_haplotype_mappings[b.first][std::get<1>(winning_pair)];
            int total_hom_kmers = barcode_hom_mappings[b.first];
            int total_kmers = kmers_per_barcode[b.first];
            int total_disagreeing_kmers = haplotype_barcode_disagree[std::get<1>(winning_pair)][b.first];
            int other = total_kmers - total_agreeing_kmers - total_hom_kmers - total_disagreeing_kmers;
            out << b.first << ": " << total_agreeing_kmers << ", " << total_disagreeing_kmers << ", " << total_hom_kmers <<", " << other << ", " << total_kmers << std::endl;
        }
    }
}


void HaplotypeScorer::add_barcode_vote(std::string barcode, std::string edge, int kmers){
    barcodes.insert(barcode);
    if (edge_haplotype_dict.find(edge) != edge_haplotype_dict.end()){
        // we only care about mappings to edges in bubbles, which will all have a key in the edge dict
        barcode_edge_mappings[barcode][edge] += kmers;
    } else if (std::find(graph.edges.begin(), graph.edges.end(), edge) != graph.edges.end()){
        barcode_hom_mappings[barcode] += kmers;
    }
    kmers_per_barcode[barcode] += kmers;

}

void HaplotypeScorer::load_mappings_from_dict(std::map<std::string, std::map<std::string, int> > & mappings) {
    for (auto edge:graph.edges){
        for (auto barcode: mappings[edge]){
            add_barcode_vote(barcode.first, edge, barcode.second);
        }
    }
}
