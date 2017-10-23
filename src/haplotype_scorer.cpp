//
// Created by Katie Barr (EI) on 15/10/2017.
//

#include <numeric>
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
/*
void print_stuff_vect_int(std::vector<int> in){
    auto a = accumulator<float,double>(); // Generate summary statistics for the offset distribution
    for (auto element; in) {
        a(element); // Call once per value
    }
    std::cout << a; // print statistics
}*/

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
    auto support_max = std::max_element(haplotype_support_vals.begin(), haplotype_support_vals.end());
    auto support_mean = avg(haplotype_support_vals);
    auto overall_support_max = std::max_element(haplotype_overall_support_vals.begin(),
                                                haplotype_overall_support_vals.end());
    auto overall_support_mean = avg(haplotype_support_vals);

    auto support_stdev = stdev(haplotype_support_vals, support_mean);

    auto overall_stdev = stdev(haplotype_overall_support_vals, overall_support_mean);
    std::cout << "Haplotype support max : " << *support_max << " min: "
              << *std::min_element(haplotype_support_vals.begin(), haplotype_support_vals.end()) << " mean: "
              << support_mean << " stdev: " << support_stdev << std::endl;
    std::cout << "Haplotype not support max : "
              << *std::max_element(haplotype_not_support_vals.begin(), haplotype_not_support_vals.end()) << " min: "
              << *std::min_element(haplotype_not_support_vals.begin(), haplotype_not_support_vals.end()) << " mean: "
              << avg(haplotype_not_support_vals) << std::endl;
    std::cout << "Haplotype overall support max: " << *overall_support_max << " min: "
              << *std::min_element(haplotype_overall_support_vals.begin(), haplotype_overall_support_vals.end())
              << " mean: " << overall_support_mean << " stdev: " << overall_stdev << std::endl;
    // max, min, mean, stdev for support, not support, overall support
    out << *support_max << "\t" << *std::min_element(haplotype_support_vals.begin(), haplotype_support_vals.end()) << "\t"
                                              << support_mean << "\t" << support_stdev << "\t" << *std::max_element(haplotype_not_support_vals.begin(), haplotype_not_support_vals.end()) << "\t"
                                                                                               << *std::min_element(haplotype_not_support_vals.begin(), haplotype_not_support_vals.end()) << "\t"
                                                                                               << avg(haplotype_not_support_vals) << "\t" << *overall_support_max << "\t"
                                                                                                                                  << *std::min_element(haplotype_overall_support_vals.begin(), haplotype_overall_support_vals.end()) << std::endl;
    for (int i=0; i < 3; i++) {
        out << "Barcode support\tHap" << "\t" << i << "\t";
        for (auto h: possible_haplotypes[supports[i].first]){
            out << h << "\t";
        }
        out << std::endl;
    }
    for (int i=0; i < 3; i++) {
        out << "Overall support\tHap" << "\t" << i << "\t";

        for (auto h: possible_haplotypes[overall_supports[i].first]){
            out << h << "\t";
        }
        out << std::endl;
    }
    max_overall_support = *support_max;
    mean_overall_support = support_mean;
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
    auto pair_support_max = std::max_element(hap_pair_support_values.begin(), hap_pair_support_values.end());
    auto not_pair_support_max = std::max_element(hap_pair_not_support_values.begin(),
                                                 hap_pair_not_support_values.end());
    auto overall_pair_support_max = std::max_element(hap_pair_support_total_score_values.begin(),
                                                     hap_pair_support_total_score_values.end());
    auto pair_support_min = std::min_element(hap_pair_support_values.begin(), hap_pair_support_values.end());
    auto not_pair_support_min = std::min_element(hap_pair_not_support_values.begin(),
                                                 hap_pair_not_support_values.end());
    auto overall_pair_support_min = std::min_element(hap_pair_support_total_score_values.begin(),
                                                     hap_pair_support_total_score_values.end());
    auto pair_support_mean = avg(hap_pair_support_values);
    auto pair_support_stdev = stdev(hap_pair_support_values, pair_support_mean);
    auto pair_not_support_mean = avg(hap_pair_not_support_values);
    auto pair_not_support_stdev = stdev(hap_pair_not_support_values, pair_not_support_mean);
    auto pair_overall_support_mean = avg(hap_pair_support_total_score_values);
    auto pair_overall_support_stdev = stdev(hap_pair_support_total_score_values, pair_overall_support_mean);
    std::cout << "pair support size: " << pair_supports.size() << " total: "
              << hap_pair_support_total_score_values.size() << std::endl;
    max_overall_pair_support = *overall_pair_support_max;
    mean_overall_pair_support = pair_overall_support_mean;
    std::cout << "Haplotype pair support max : " << *pair_support_max << "min : " << *pair_support_min << " mean: "
              << pair_support_mean << " stdev: " << pair_support_stdev << std::endl;
    std::cout << "Haplotype pair not support max : " << *not_pair_support_max << "min : " << *not_pair_support_min
              << " mean: " << pair_not_support_mean << " stdev: " << pair_not_support_stdev << std::endl;
    std::cout << "Haplotype pair overall support max : " << *overall_pair_support_max << "min : "
              << *overall_pair_support_min << " mean: " << pair_overall_support_mean << " stdev: "
              << pair_overall_support_stdev << std::endl;
    // max, min, mean, stdev for support, not support, overall support
    out << *pair_support_max << "\t" << *std::min_element(hap_pair_support_values.begin(), hap_pair_support_values.end()) << "\t"
            << avg(hap_pair_support_values) << "\t" << pair_support_stdev << "\t" << *std::max_element(hap_pair_not_support_values.begin(), hap_pair_not_support_values.end()) << "\t"
                                                     << *std::min_element(hap_pair_not_support_values.begin(), hap_pair_not_support_values.end()) << "\t"
                                                     << avg(hap_pair_not_support_values) << "\t" << stdev(hap_pair_not_support_values, avg(hap_pair_not_support_values))<< "\t"<< *overall_pair_support_max << "\t"
            << *std::min_element(hap_pair_support_total_score_values.begin(), hap_pair_support_total_score_values.end()) << std::endl;

    for (int i=0; i < 3; i++) {
        out << "Barcode support\tHapP1" << "\t" << i << "\t";

        for (auto h: possible_haplotypes[std::get<0>(pair_supports[i].first)]){
            out << h << "\t";
        }
        out << std::endl;
    }
    for (int i=0; i < 3; i++) {
        out << "Barcode support\tHapP2" << "\t" << i << "\t";
        for (auto h: possible_haplotypes[std::get<1>(pair_supports[i].first)]){
            out << h << "\t";
        }
        out << std::endl;
    }
    for (int i=0; i < 3; i++) {
        out << "Overall support\tHapP1" << "\t" << i << "\t";

        for (auto h: possible_haplotypes[std::get<0>(pair_overall_supports[i].first)]){
            out << h << "\t";
        }
        out << std::endl;
    }
    for (int i=0; i < 3; i++) {
        out << "Barcode support\tHapP2" << "\t" << i << "\t";

        for (auto h: possible_haplotypes[std::get<1>(pair_overall_supports[i].first)]){
            out << h << "\t";
        }
        out << std::endl;
    }
}

int HaplotypeScorer::score_haplotypes(std::string outfile) {
    //initialize score arrays- index is haplotype index
    int haplotype_support[possible_haplotypes.size()] = {0};
    int haplotype_not_support[possible_haplotypes.size()] = {0};
    int haplotype_overall_support[possible_haplotypes.size()] = {0};
    std::map<std::pair<int, int>, int> hap_pair_not_support;
    std::map<std::pair<int, int>, int> hap_pair_support;
    std::map<std::pair<int, int>, int> hap_pair_support_total_score;
    std::string barcode;
    std::ofstream mappings("barcode_mappings.txt");
    for (auto &bm: barcode_haplotype_mappings) {
        for (auto e:bm.second) {
            mappings << bm.first << " : " << e.first << " " << e.second << std::endl;
        }
    }
    for (auto &bm: barcode_haplotype_mappings) {
        barcode = bm.first;
        // winner_for_this_barcode = [h for h in self.barcode_mappings[barcode] if self.barcode_mappings[barcode][h] == np.max(hap_support_dict.values())]
        std::vector<int> winners = winner_for_barcode(barcode); // ideally should be length 1
        // for haplotype in range(len(self.list_of_possible_haplotypes)/2):
        for (int hap = 0; hap < possible_haplotypes.size() / 2; hap++) {
            // pair = len(self.list_of_possible_haplotypes) -1 -haplotype
            int pair = possible_haplotypes.size() - 1 - hap;
            if (bm.second.find(hap) != bm.second.end()) {
                haplotype_overall_support[hap] += bm.second[hap];
                hap_pair_support_total_score[std::make_pair(hap, pair)] += bm.second[hap];
                if (winners.size() > 0) {
                    if (std::find(winners.begin(), winners.end(), hap) != winners.end()) {
                        haplotype_support[hap] += 1;
                        hap_pair_support[std::make_pair(hap, pair)] += 1;
                        haplotype_barcode_agree[hap][barcode] += bm.second[hap];
                        haplotype_barcode_disagree[hap][barcode] += bm.second[pair];
                    }
                }
            }

            if (bm.second.find(pair) != bm.second.end()) {
                haplotype_overall_support[pair] += bm.second[pair];
                hap_pair_support_total_score[std::make_pair(hap, pair)] += bm.second[pair];
                if (winners.size() > 0) {
                    if (std::find(winners.begin(), winners.end(), pair) != winners.end()) {
                        haplotype_support[pair] += 1;
                        hap_pair_support[std::make_pair(hap, pair)] += 1;
                        // for each barcode which selects the winner, need total kmers agreeing/disagreeing
                        haplotype_barcode_agree[pair][barcode] += bm.second[pair];
                        // for each barcode, and each candidate haplotype, need to find out how many kmers from that don't support that haplotype- which is number of kmers for that barcode that support other in pair
                        haplotype_barcode_disagree[pair][barcode] += bm.second[hap];
                    }
                }
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
    std::cout << "Calculating barcode haplotype support for " << barcode_edge_mappings.size() << " mappings"<< std::endl;
    for (auto &mapping:barcode_edge_mappings){
        //std::cout << "Checking barcode " << mapping.first <<std::endl;
        // if barcode maps to more than 1 edge in bubbles and maximum support is greater than 1
        //auto edge_support_max = std::max_element(std::begin(mapping.second), std::end(mapping.second), [] ( std::map<std::string, int> &p1,  std::map<std::string, int> &p2) {return p1.second < p2.second});
        // if len(self.barcode_edge_mappings[barcode].keys()) > 1:
        if (mapping.second.size() > 1){
            //edge_support = {edge:self.barcode_edge_mappings[barcode][edge] for edge in self.barcode_edge_mappings[barcode] if edge in self.graph.edge_bubble_dict.keys()}
            std::vector<std::string> edges;
            std::vector<int> scores;
            for (auto e: mapping.second){
                edges.push_back(e.first);
                scores.push_back(e.second);
            }
            // elif len(edges) > 1 and np.max(edge_support.values()) > 1:
            if (*std::max_element(scores.begin(), scores.end())> 1) {
                //  for i, haplotype in enumerate(self.list_of_possible_haplotypes)
                for (int i = 0; i < possible_haplotypes.size(); i++) {
                    std::vector<std::string> edges_in_haplotype;
                    std::vector<std::string> h;
                    h = possible_haplotypes[i];
                    // find all edges in each haplotype that this barcode maps to
                    for (auto e1: edges) {
                        //edges_in_haplotype = [e for e in haplotype if e in edges]
                        if (std::find(h.begin(), h.end(), e1) != h.end()) {
                            edges_in_haplotype.push_back(e1);
                        }
                    }
                    // somewhat arbitrary rule to decide if the barcode supports a haplotype enough
                    // if len(edges_in_haplotype)>= len(edges)/2 and len(edges_in_haplotype) > 1:
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
        //std::cout << "barcode " << mapping.first << " supports " << haplotypes_supported << std::endl;

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
