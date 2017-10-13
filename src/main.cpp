#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

#include <sys/stat.h>

//#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

#include "../deps/cxxopts/include/cxxopts.hpp"
#include "graph.h"

Graph load_subgraph(std::string subgraph_filename){
    ifstream subgraph(subgraph_filename);

}

bool check_or_create_directory(std::string &output_prefix) {
    if (output_prefix.back() != '/') {
        output_prefix.push_back('/');
    }
    struct stat sb{};
    bool validate_dir(false);
    if (stat(output_prefix.c_str(), &sb) != 0) {
        if (errno == ENOENT) {
            mode_t mask = umask(0);
            umask(mask);
            mkdir(output_prefix.c_str(), mode_t(0777 - mask));
            validate_dir = true;
        }
        if (stat(output_prefix.c_str(), &sb) != 0) {
            perror(output_prefix.c_str());
            validate_dir = false;
        }
    } else if (!S_ISDIR(sb.st_mode)) {
        std::cout << output_prefix << " is not a directory " << std::endl;
    } else {
        validate_dir = true;
    }
    return validate_dir;
}

int main(int argc, char **argv) {

    const unsigned int GB(1024*1024*1024);
    unsigned int mem_limit(10);

    std::string asm_filename;
    std::string fastq_filename;
    uint16_t min_count(0);
    uint32_t max_count(100000);
    uint16_t max_coverage(1);
    uint32_t min_read_length(1000), min_contig_length(1000);
    uint32_t min_kmers_to_call_match(10);
    uint32_t min_seen_contig_to_write_output(1);

    std::string set_filelist;
    std::string output_prefix;

    uint8_t m(1),n(1),k(31);

    try
    {
        std::time_t t = std::time(nullptr);
        std::tm tm = *std::localtime(&t);
        std::string outdefault(
                std::to_string(tm.tm_year + 1900) + '-' + std::to_string(tm.tm_mon) + '-' + std::to_string(tm.tm_mday) +
                '_' + std::to_string(tm.tm_hour) + std::to_string(tm.tm_min));

        cxxopts::Options options("seq-sorter", "Sequence linking tool using long/linked reads.");

        options.add_options()("help", "Print help")("a,assembly", "Sequence reference to link",
                                                    cxxopts::value<std::string>(asm_filename), "FASTA - Sequence file")(
                "r,long_reads", "Reads to generate sequence-to-sequence links",
                cxxopts::value<std::string>(fastq_filename), "FASTQ - Reads")("min_count",
                                                                              "Minimum count to consider a link",
                                                                              cxxopts::value<uint16_t>(
                                                                                      min_count)->default_value("0"),
                                                                              "uint")("max_count",
                                                                                      "Maximum count to consider a link",
                                                                                      cxxopts::value<uint32_t>(
                                                                                              max_count)->default_value(
                                                                                              "100000"), "uint")(
                "max_coverage", "max coverage for a kmer to be considered",
                cxxopts::value<uint16_t>(max_coverage)->default_value("1"), "uint")("min_read_length",
                                                                                    "minimum contig length",
                                                                                    cxxopts::value<uint32_t>(
                                                                                            min_read_length)->default_value(
                                                                                            "1000"), "uint")(
                "min_contig_length", "minimum read length",
                cxxopts::value<uint32_t>(min_contig_length)->default_value("1000"), "uint")("min_kmers_to_call_match",
                                                                                            "minimum number of kmers to call a Read->Contig match",
                                                                                            cxxopts::value<uint32_t>(
                                                                                                    min_kmers_to_call_match)->default_value(
                                                                                                    "10"), "uint")(
                "o,output", "output file prefix", cxxopts::value<std::string>(output_prefix)->default_value(outdefault),
                "prefix_dir");

        options.add_options("Output options")("min_seen_contig_to_write_output",
                                              "minimum number of seen contigs to report read on output (1)",
                                              cxxopts::value<uint32_t>(min_seen_contig_to_write_output));

        options.add_options("Skip-mer shape (m every n, total k)")
                ("m,used_bases", "m (1)", cxxopts::value<uint8_t>(m))
                ("n,skipped_bases", "n (1)", cxxopts::value<uint8_t>(n))
                ("k,total_bases", "k (31)", cxxopts::value<uint8_t>(k));

        options.add_options("Performance")("mem_limit", "Memory limit in GB",
                                           cxxopts::value<unsigned int>(mem_limit)->default_value("10"));

        options.parse(argc, argv);

        if (0 != options.count("help"))
        {
            std::cout << options.help({"", "Performance"}) << std::endl;
            exit(0);
        }
}        if (options.count("o")!=1 /*or options.count("i")<2*/) {
        std::cout << "Error: please specify input files and output prefix"<<std::endl
                  <<" Use option --help to check command line arguments." << std::endl;
        exit(1);
    }


    if (asm_filename.empty()) {
        std::cout << "Error: The assembly file parameter wasn't specified, " << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
    if (fastq_filename.empty()) {
        std::cout << "Error: The read file parameter wasn't specified, " << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
} catch (const cxxopts::OptionException& e)
{
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit(1);
}