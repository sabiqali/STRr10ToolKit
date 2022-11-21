#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <stdbool.h>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <unordered_set>
#include <chrono>
#include <sstream>
#include <getopt.h>
#include "htslib/sam.h"
#include "htslib/bam.h"
#include "./str_utilities.h"

KSEQ_INIT(gzFile, gzread)

static const char *ALIGN_MESSAGE = 
"Usage: ./strr10toolkit [OPTIONS] --bam input.bam --reference reference_genome.fasta --output_file_name output_file_name --output_directory output_directory\n"
"Align reads to control oligos provided.\n"
"\n"
"  -b, --bam=FILE                       the input bam file\n"
"  -r, --reference=FILE                 the reference genome\n"
"  -o, --output_file_name=FILE          the output file name without extension\n"
"  -d, --output_directory=PATH          the directory where all output files will be generated\n"
"      --min_ins_size=NUM               the minimum number of bases in the insert to call as an STR. Default Minimum is 50bp.\n"
"      --is_phased                      is phasing data present in the bam?\n"
"      --min_read_support=NUM           the minimum number of reads that support the insert call. Minimum is set at 3 by default.\n"
"      --chromosomes=FILE               the chromosomes to be searched for STRs. Alternatively, you can provide a list of chromosomes in text file separated by a whitespace\n"
"  -c, --clean                          cleanup any intermediate files\n"
"      --window_size=NUM                the size of the search window. Defaults to 5000bp\n"
"      --min_repeat_size=NUM            minimum length of repeat motif. Defaults to 3\n"
"      --max_repeat_size=NUM            maximum length of repeat motif. Defaults to 6\n"
"      --min_map_quality=NUM            minimum mapping quality of individual reads. Defaults to 20\n"
"      --discovery_sensitivity=NU<      minimum number of repeat units in the repeat expansion to be considered, to be used with min_ins_size. Defaults to 10\n"
"  -t, --threads=NUM                    use NUM threads. Defaults to 1\n";

namespace opt {
    static std::string bam_file;
    static std::string reference_file;
    static std::string output_file_name;
    static std::string output_directory;
    static std::string chromosomes;
    static int min_ins_size = 50;
    static int is_phased = 0;
    static int min_read_support = 3;
    static int clean = 0;
    static int num_threads = 1;
    static int window_size = 5000;
    static int min_repeat_size = 3;
    static int max_repeat_size = 6;
    static int min_map_quality = 20;
    static int discovery_sensitivity = 20;
}

static const char* shortopts = "b:r:o:d:c:t";

enum { OPT_HELP = 1, OPT_VERSION, OPT_CHROMOSOME , OPT_MIN_INS_SIZE, OPT_IS_PHASED, OPT_MIN_READ_SUPPORT, OPT_WINDOW_SIZE, OPT_MIN_REPEAT_SIZE, OPT_MAX_REPEAT_SIZE, OPT_MIN_MAP_QUALITY, OPT_DISCOVERY_SENSITIVITY };

static const struct option longopts[] = {
    { "bam_file",             required_argument, NULL, 'b' },
    { "reference_file",       required_argument, NULL, 'r' },
    { "output_file_name",     required_argument, NULL, 'o' },
    { "output_directory",     required_argument, NULL, 'd' },
    { "chromosomes",          required_argument, NULL, OPT_CHROMOSOME },
    { "min_ins_size",         required_argument, NULL, OPT_MIN_INS_SIZE },
    { "is_phased",            no_argument,       NULL, OPT_IS_PHASED },
    { "min_read_support",     required_argument, NULL, OPT_MIN_READ_SUPPORT },
    { "clean",                no_argument,       NULL, 'c' },
    { "num_threads",          required_argument, NULL, 't' },
    { "window_size",          required_argument, NULL, OPT_WINDOW_SIZE },
    { "min_repeat_size",      required_argument, NULL, OPT_MIN_REPEAT_SIZE },
    { "max_repeat_size",      required_argument, NULL, OPT_MAX_REPEAT_SIZE },
    { "min_map_quality",      required_argument, NULL, OPT_MIN_MAP_QUALITY },
    { "discovery_sensitivity",required_argument, NULL, OPT_DISCOVERY_SENSITIVITY },
    { "help",                 no_argument,       NULL, OPT_HELP },
    { "version",              no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

struct per_read_struct {
    int region_start;
    int region_end;
    int region_ref_start;
    int region_ref_end;
    std::string motif;
    int haplotype;
    int size;
    float min_methylation;
    float max_methylation;
    float avg_methylation;
    std::string interruption_motif;
    char strand;
    std::string query_name;
};

struct per_window_struct {
    std::vector<per_read_struct> window_aggregate;
    std::map<std::string,int> motif_aggregate;
};

void parse_align_options(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'b': arg >> opt::bam_file; break;
            case 'r': arg >> opt::reference_file; break;
            case 'o': arg >> opt::output_file_name; break;
            case 'd': arg >> opt::output_directory; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'c': opt::clean = 1; break;
            case OPT_CHROMOSOME: arg >> opt::chromosomes; break;
            case OPT_MIN_INS_SIZE: arg >> opt::min_ins_size; break;
            case OPT_IS_PHASED: opt::is_phased = 1; break;
            case OPT_MIN_READ_SUPPORT: arg >> opt::min_read_support; break;
            case OPT_WINDOW_SIZE: arg >> opt::window_size; break;
            case OPT_MIN_REPEAT_SIZE: arg >> opt::min_repeat_size; break;
            case OPT_MAX_REPEAT_SIZE: arg >> opt::max_repeat_size; break;
            case OPT_MIN_MAP_QUALITY: arg >> opt::min_map_quality; break;
            case OPT_DISCOVERY_SENSITIVITY: arg >> opt::discovery_sensitivity; break;
            case OPT_HELP:
                std::cout << ALIGN_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << ALIGN_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if(opt::num_threads <= 0) {
        std::cerr << "STRr10ToolKit : invalid number of threads: " << opt::num_threads << "\n";
        die = true;
    }

    if(opt::bam_file.empty()) {
        std::cerr << "STRr10ToolKit : an input bam file must be provided\n";
        die = true;
    }

    if(opt::reference_file.empty()) {
        std::cerr << "STRr10ToolKit : a reference must be provided\n";
        die = true;
    }

    if(opt::chromosomes.empty()) {
        std::cerr << "STRr10ToolKit : the chromosomes must be specified\n";
        die = true;
    }

    if(opt::output_file_name.empty() || opt::output_directory.empty()) {
        std::cerr << "STRr10ToolKit : an output file and directory must be provided\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << ALIGN_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int ref_chr_size(std::string chr_name) {
    switch (chr_name)
    {
    case "chr1":
        return 248956422;
    case "chr2":
        return 242193529;
    case "chr3":
        return 198295559;
    case "chr4":
        return 190214555;
    case "chr5":
        return 181538259;
    case "chr6":
        return 170805979;
    case "chr7":
        return 159345973;
    case "chr8":
        return 145138636;
    case "chr9":
        return 138394717;
    case "chr10":
        return 133797422;
    case "chr11":
        return 135086622;
    case "chr12":
        return 133275309;
    case "chr13":
        return 114364328;
    case "chr14":
        return 107043718;
    case "chr15":
        return 101991189;
    case "chr16":
        return 90338345;
    case "chr17":
        return 83257441;
    case "chr18":
        return 80373285;
    case "chr19":
        return 58617616;
    case "chr20":
        return 64444167;
    case "chr21":
        return 46709983;
    case "chr22":
        return 50818468;
    case "chrX":
        return 156040895;
    case "chrY":
        return 57227415;
    default:
        return 0;
    }
}

int main(int argc, char *argv[])  {

    parse_align_options(argc , argv);

    //TODO::Init the bam file and get reads in window in chromosomes requested. Pass the alignments to the utilities to compute the stats
    samFile *fp = sam_open(opt::bam_file, "r");
    std::string bai_file = opt::bam_file + ".bai";
    hts_idx_t *idx = sam_index_load(fp, bai_file);
    bam_hdr_t *h = sam_hdr_read(fp);
    bam1_t *b = bam_init1();

    //get the chr list from file
    std::ifstream chr_file(opt::chromosomes);
    std::string chr_list_line; 
    std::string token;
    std::vector<std::string> chr_list;
    std::getline(chr_file, chr_list_line); //only 1 line in file is expected

    std::stringstream tmp_chr_list(chr_list_line);
    while (getline(tmp_chr_list, token, ' ')){
        chr_list.push_back(token);
    }

    for (auto chr: chr_list) {
        int lower_limit = 0;
        int upper_limit = opt::window_size;

        while(upper_limit <= ref_chr_size(chr)) {
            per_window_struct* window_output = new per_read_struct();

            std::string region = chr+":"+std::to_string(lower_limit)+"-"+std::to_string(upper_limit);
            //generates iterator over region
            hts_itr_t *itr = sam_itr_querys(idx, h, region);
            /* or do i use this?:
            int bam_fetch(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func);
            */
            while (sam_itr_next(fp, itr, b) >= 0) {
                per_read_struct* read_output = new per_read_struct();

                //manip reads in region and send to the appropriate functions to get the results
                if (b->core.tid < 0) continue;

                int read_pos_counter = 0;
                int ref_pos_counter = 0;
                std::string query_sequence = bam1_seq(b);

                auto ref_start_pos = b->core.pos;

                /*printf("%s\t%d\t%d\n", h->target_name[b->core.tid], b->core.pos,
                b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)));*/
                //bam handling derived from: https://www.biostars.org/p/4211/
                int depth[REF_LEN], x, j, k;
                uint32_t *cigar = bam1_cigar(b);
                for (k = 0, x = b->core.pos; k < b->core.n_cigar; ++k) {
                    int op = cigar[k]&16;
                    int l = cigar[k]>>4;
                    switch(op) {
                        case BAM_CMATCH:
                            //printf("M");
                            read_pos_counter += l;
                            ref_pos_counter += l;
                            break;
                        case BAM_CHARD_CLIP:
                            //printf("H");
                            break;
                        case BAM_CSOFT_CLIP:
                            //printf("S");
                            if(l<opt::min_ins_size){
                                read_pos_counter += l;
                            }
                            else {
                                //call functions here
                                std::string sequence_of_interest = query_sequence[read_pos_counter:l];

                                auto [potential_str_sequence,potential_count_from_discovery] = decompose_string(sequence_of_interest,opt::min_repeat_size,opt::max_repeat_size);

                                auto [count,interrupting_motif] = detect_size(sequence_of_interest,potential_str_sequence);

                                auto [min_methylation,max_methylation,avg_methylation] = detect_methylation(read_pos_counter,l,b); 

                                if(count<=opt::discovery_sensitivity) {
                                    read_pos_counter += l;
                                    break;
                                }
                                if(potential_str_sequence.empty()) {
                                    read_pos_counter += l;
                                    break;
                                }

                                read_output->region_ref_start = ref_start_pos + ref_pos_counter;
                                read_output->region_ref_end = read_output->region_ref_start + 1;
                                read_output->region_start = read_pos_counter;
                                read_output->region_end = read_pos_counter+l;
                                read_output->interruption_motif = interrupting_motif;
                                read_output->size = count;
                                read_output->query_name = bam1_qname(b);
                                if(bam_is_rev(b)) {
                                    read_output->motif = dna_reverse_complement(potential_str_sequence);
                                    read_output->strand = '-';
                                } 
                                else {
                                    read_output->motif = potential_str_sequence;
                                    read_output->strand = '+';
                                }
                                read_output->min_methylation = min_methylation;
                                read_output->max_methylation = max_methylation;
                                read_output->avg_methylation = avg_methylation;

                                read_pos_counter += l;
                            }
                            break;
                        case BAM_CDEL:
                            //printf("D"); does this span positions on the read???
                            ref_pos_counter += l;
                            break;
                        case BAM_CPAD:
                            //printf("P");
                            read_pos_counter += l;
                            break;
                        case BAM_CINS:
                            //printf("I");
                            if(l<opt::min_ins_size){
                                read_pos_counter += l;
                            }
                            else {
                                //call functions here
                                std::string sequence_of_interest = query_sequence[read_pos_counter:l];

                                auto [potential_str_sequence,potential_count_from_discovery] = decompose_string(sequence_of_interest,opt::min_repeat_size,opt::max_repeat_size);

                                auto [count,interrupting_motif] = detect_size(sequence_of_interest,potential_str_sequence);

                                auto [min_methylation,max_methylation,avg_methylation] = detect_methylation(read_pos_counter,l,b); 

                                if(count<=opt::discovery_sensitivity) {
                                    read_pos_counter += l;
                                    break;
                                }
                                if(potential_str_sequence.empty()) {
                                    read_pos_counter += l;
                                    break;
                                }

                                //auto ref_region_start_pos = bam_calend(b->core,cigar);
                                //auto read_size = b->core.l_qseq;

                                read_output->region_ref_start = ref_start_pos + ref_pos_counter;
                                read_output->region_ref_end = read_output->region_ref_start + 1;
                                read_output->region_start = read_pos_counter;
                                read_output->region_end = read_pos_counter+l;
                                read_output->interruption_motif = interrupting_motif;
                                read_output->size = count;
                                read_output->query_name = bam1_qname(b);
                                if(bam_is_rev(b)) {
                                    read_output->motif = dna_reverse_complement(potential_str_sequence);
                                    read_output->strand = '-';
                                } 
                                else {
                                    read_output->motif = potential_str_sequence;
                                    read_output->strand = '+';
                                }
                                read_output->min_methylation = min_methylation;
                                read_output->max_methylation = max_methylation;
                                read_output->avg_methylation = avg_methylation;

                                read_pos_counter += l;
                            }
                            break;
                        case BAM_CREF_SKIP:
                            //printf("S");
                            if(l<opt::min_ins_size){
                                read_pos_counter += l;
                            }
                            else {
                                //call functions here
                                std::string sequence_of_interest = query_sequence[read_pos_counter:l];

                                auto [potential_str_sequence,potential_count_from_discovery] = decompose_string(sequence_of_interest,opt::min_repeat_size,opt::max_repeat_size);

                                auto [count,interrupting_motif] = detect_size(sequence_of_interest,potential_str_sequence);

                                auto [min_methylation,max_methylation,avg_methylation] = detect_methylation(read_pos_counter,l,b); 

                                if(count<=opt::discovery_sensitivity) {
                                    read_pos_counter += l;
                                    break;
                                }
                                if(potential_str_sequence.empty()) {
                                    read_pos_counter += l;
                                    break;
                                }

                                read_output->region_ref_start = ref_start_pos + ref_pos_counter;
                                read_output->region_ref_end = read_output->region_ref_start + 1;
                                read_output->region_start = read_pos_counter;
                                read_output->region_end = read_pos_counter+l;
                                read_output->interruption_motif = interrupting_motif;
                                read_output->size = count;
                                read_output->query_name = bam1_qname(b);
                                if(bam_is_rev(b)) {
                                    read_output->motif = dna_reverse_complement(potential_str_sequence);
                                    read_output->strand = '-';
                                } 
                                else {
                                    read_output->motif = potential_str_sequence;
                                    read_output->strand = '+';
                                }
                                read_output->min_methylation = min_methylation;
                                read_output->max_methylation = max_methylation;
                                read_output->avg_methylation = avg_methylation;

                                read_pos_counter += l;
                            }
                            break;
                        default:
                            //printf("?");
                            std::cerr<<"STRr10ToolKit::Cigar_Parse: cannot parse cigar element";
                    }
                }
                window_output->window_aggregate.push_back(read_output);
                if(window_output->motif_aggregate.find(read_output->motif) == window_output->motif_aggregate.end()) {
                    window_output->motif_aggregate.insert({read_output->motif,1});
                }
                else {
                    window_output->motif_aggregate[read_output->motif] += 1;
                }
                int max_read_support=0;
                for(auto &entry: window_output->motif_aggregate) {
                    if(entry->second > max_read_support) {
                        max_read_support = entry->second;
                    }
                }
                if(max_read_support >= opt::min_read_support) {
                    //print out the stats from this window
                }
                //otherwise, there aren't any STRs that pass all the filters. moving to the next window
            }
            hts_itr_destroy(itr);
        }
    }

    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fp);
}