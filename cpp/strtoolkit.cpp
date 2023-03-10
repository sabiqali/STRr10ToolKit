#include <zlib.h>
#include <stdio.h>
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
//#include "./htslib/bam.h"
#include "htslib/htslib/vcf.h"
#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"
#include "htslib/htslib/hts.h"
#include "str_utilities.h"
#include "bloom_filter.h"

static const char *TOOLKIT_MESSAGE = 
"Usage: ./strr10toolkit [OPTIONS] --bam input.bam --reference reference_genome.fasta --output_file_name output_file_name.tsv --output_directory /path/to/output/directory\n"
"Toolkit to detect and analyse STR loci\n"
"\n"
"  -b, --bam=FILE                       the input bam file\n"
"  -r, --reference=FILE                 the reference genome\n"
"  -o, --output_file_name=FILE          the output file name without extension\n"
"  -d, --output_directory=PATH          the directory where all output files will be generated\n"
"      --min_ins_size=NUM               the minimum number of bases in the insert event to call as an STR.(Default: 50bp)\n"
"      --is_phased                      is phasing data present in the bam?\n"
"  -v, --verbose                        will output the data from intermediate steps. Requires output_file_name and output_directory to be set, if used. will print to stdout if not set\n"
"      --min_read_support=NUM           the minimum number of reads that support the insert call.(Default: 3)\n"
"      --chromosomes=FILE               the chromosomes to be searched for STRs. Alternatively, you can provide a list of chromosomes in text file separated by a whitespace\n"
"  -c, --clean                          cleanup any intermediate files\n"
"      --window_size=NUM                the size of the search window.(Default: 5000bp)\n"
"      --min_repeat_size=NUM            minimum length of repeat motif.(Default: 3)\n"
"      --max_repeat_size=NUM            maximum length of repeat motif. (Default: 6)\n"
"      --min_map_quality=NUM            minimum mapping quality of individual reads.(Default: 20)\n"
"      --discovery_sensitivity=NUM      minimum number of repeat units in the repeat expansion to be considered, to be used with min_ins_size.(Default: 10)\n"
"  -t, --threads=NUM                    use NUM threads.(Default: 1)\n";

namespace opt {
    static std::string bam_file;
    static std::string ref_file;
    static std::string output_file;
    static std::string output_directory;
    static std::string chromosome_file;
    static int min_ins_size = 50;
    static int is_phased = 0;
    static int verbose = 0;
    static int min_read_support = 3;
    static int clean_flag = 0 ;
    static int window_size = 5000;
    static int min_repeat_size = 3;
    static int max_repeat_size = 6;
    static int min_map_quality = 20;
    static int discovery_sensitivity = 10;
    static int num_threads = 1;
}

static const char* shortopts = "b:r:o:d:v:c:t";

enum { OPT_HELP = 1, OPT_VERSION, OPT_CHROMOSOME_FILE, OPT_MIN_INS_SIZE, OPT_IS_PHASED, OPT_MIN_READ_SUPPORT, OPT_WINDOW_SIZE, OPT_MIN_REPEAT_SIZE, OPT_MAX_REPEAT_SIZE, OPT_MIN_MAP_QUALITY, OPT_DISCOVERY_SENSITIVITY };

static const struct option longopts[] = {
    { "bam_file",              required_argument, NULL, 'b' },
    { "ref_file",              required_argument, NULL, 'r' },
    { "output_file",           required_argument, NULL, 'o' },
    { "output_directory",      required_argument, NULL, 'd' },
    { "chromosome",            required_argument, NULL, OPT_CHROMOSOME_FILE },
    { "min_ins_size",          required_argument, NULL, OPT_MIN_INS_SIZE },
    { "is_phased",             no_argument,       NULL, OPT_IS_PHASED },
    { "verbose",               no_argument,       NULL, 'v' },
    { "min_read_support",      required_argument, NULL, OPT_MIN_READ_SUPPORT },
    { "clean",                 no_argument,       NULL, 'c' },
    { "window_size",           required_argument, NULL, OPT_WINDOW_SIZE },
    { "min_repeat_size",       required_argument, NULL, OPT_MIN_REPEAT_SIZE },
    { "max_repeat_size",       required_argument, NULL, OPT_MAX_REPEAT_SIZE },
    { "min_map_quality",       required_argument, NULL, OPT_MIN_MAP_QUALITY },
    { "discovery_sensitivity", required_argument, NULL, OPT_DISCOVERY_SENSITIVITY },
    { "threads",               required_argument, NULL, 't' },
    { "help",                  no_argument,       NULL, OPT_HELP },
    { "version",               no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

struct per_read_struct {
    uint32_t region_start;
    uint32_t region_end;
    uint32_t region_ref_start;
    uint32_t region_ref_end;
    std::string motif;
    int haplotype;
    int size;
    float min_methylation;
    float max_methylation;
    float avg_methylation;
    float up_min_methylation;
    float up_max_methylation;
    float up_avg_methylation;
    float down_min_methylation;
    float down_max_methylation;
    float down_avg_methylation;
    std::string interruption_motif;
    char strand;
    std::string query_name;
};

struct per_window_struct {
    std::vector<per_read_struct> window_aggregate;
    std::map<std::string,int> motif_aggregate;
};

void parse_align_options2(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'b': arg >> opt::bam_file; break;
            case 'r': arg >> opt::ref_file; break;
            case 'o': arg >> opt::output_file; break;
            case 'd': arg >> opt::output_directory; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'c': opt::clean_flag = 1; break;
	    case 'v': opt::verbose = 1; break;
            case OPT_CHROMOSOME_FILE: arg >> opt::chromosome_file; break;
            case OPT_MIN_INS_SIZE: arg >> opt::min_ins_size; break;
            case OPT_IS_PHASED: opt::is_phased = 1; break;
            case OPT_MIN_READ_SUPPORT: arg >> opt::min_read_support; break;
            case OPT_WINDOW_SIZE: arg >> opt::window_size; break;
            case OPT_MIN_REPEAT_SIZE: arg >> opt::min_repeat_size; break;
            case OPT_MAX_REPEAT_SIZE: arg >> opt::max_repeat_size; break;
            case OPT_MIN_MAP_QUALITY: arg >> opt::min_map_quality; break;
            case OPT_DISCOVERY_SENSITIVITY: arg >> opt::discovery_sensitivity; break;
            case OPT_HELP:
                std::cout << TOOLKIT_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << TOOLKIT_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    /*if(argc - optind > 0) {
        opt::region = argv[optind++];
    }
    if (argc - optind > 0) {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }*/

    if(opt::num_threads <= 0) {
        std::cerr << "strtoolkit : invalid number of threads: " << opt::num_threads << "\n";
        die = true;
    }

    if(opt::bam_file.empty()) {
        std::cerr << "strtoolkit: a --bam file must be provided\n";
        die = true;
    }

    if(opt::ref_file.empty()) {
        std::cerr << "strtoolkit: a --ref file must be provided\n";
        die = true;
    }

    if(opt::verbose && opt::output_file.empty() && opt::output_directory.empty()) {
        std::cerr << "strtoolkit: the output file name and the directory has to be provided for a verbose output\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << TOOLKIT_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int ref_chr_size(std::string chr_name) {

    if(chr_name.compare("chr1") == 0)
        return 248956422;
    else if(chr_name.compare("chr2") == 0)
        return 242193529;
    else if(chr_name.compare("chr3") == 0)
        return 198295559;
    else if(chr_name.compare("chr4") == 0)
        return 190214555;
    else if(chr_name.compare("chr5") == 0)
        return 181538259;
    else if(chr_name.compare("chr6") == 0)
        return 170805979;
    else if(chr_name.compare("chr7") == 0)
        return 159345973;
    else if(chr_name.compare("chr8") == 0)
        return 145138636;
    else if(chr_name.compare("chr9") == 0)
        return 138394717;
    else if(chr_name.compare("chr10") == 0)
        return 133797422;
    else if(chr_name.compare("chr11") == 0)
        return 135086622;
    else if(chr_name.compare("chr12") == 0)
        return 133275309;
    else if(chr_name.compare("chr13") == 0)
        return 114364328;
    else if(chr_name.compare("chr14") == 0)
        return 107043718;
    else if(chr_name.compare("chr15") == 0)
        return 101991189;
    else if(chr_name.compare("chr16") == 0)
        return 90338345;
    else if(chr_name.compare("chr17") == 0)
        return 83257441;
    else if(chr_name.compare("chr18") == 0)
        return 80373285;
    else if(chr_name.compare("chr19") == 0)
        return 58617616;
    else if(chr_name.compare("chr20") == 0)
        return 64444167;
    else if(chr_name.compare("chr21") == 0)
        return 46709983;
    else if(chr_name.compare("chr22") == 0)
        return 50818468;
    else if(chr_name.compare("chrX") == 0)
        return 156040895;
    else if(chr_name.compare("chrY") == 0)
        return 57227415;
    else
        return 0;
}

int main(int argc, char *argv[])  {

    parse_align_options2(argc , argv);

    //TODO::Init the bam file and get reads in window in chromosomes requested. Pass the alignments to the utilities to compute the stats
    samFile *fp = sam_open(opt::bam_file.c_str(), "r");
    std::string bai_file = opt::bam_file + ".bai";
    hts_idx_t *idx = sam_index_load(fp, bai_file.c_str());
    bam_hdr_t *h = sam_hdr_read(fp);
    bam1_t *b = bam_init1();

    //get the chr list from file
    std::ifstream chr_file(opt::chromosome_file);
    std::string chr_list_line; 
    std::string token;
    std::vector<std::string> chr_list;
    std::getline(chr_file, chr_list_line); //only 1 line in file is expected

    //init bloom filter
    bloom_filter bf_obj(220556452, 8, 20000000);
    bf_obj.allocate_data();

    std::stringstream tmp_chr_list(chr_list_line);
    while (getline(tmp_chr_list, token, ' ')){
        chr_list.push_back(token);
    }

    std::cout<<"read_name\tchromosome\tref_start\tref_end\tread_start\tread_end\tmotif\tinteruption_motif\tcount\thaplotype\tavg_methylation\tmin_methylation\tmax_methylation\tupstream_avg_methylation\tupstream_min_methylation\tupstream_max_methylation\tdownstream_avg_methylation\tdownstream_min_methylation\tdownstream_max_methylation\n";

    for (auto chr: chr_list) {
        int lower_limit = 0;
        int upper_limit = opt::window_size;

        //TODO::potentially implement a bloom filter to keep track of all reads which have been counted or visited. if you see a read which has been visited, skip it. read name is 36 chars.
        //TODO::or remove the windows and go through read by read and group by where the repeat has been found. 
        //TODO::If we go with a windowed approach, possibly label the loci ranges to differentiate the output. 
        while(upper_limit <= ref_chr_size(chr)) {
            per_window_struct* window_output = new per_window_struct();

            std::string region = chr+":"+std::to_string(lower_limit)+"-"+std::to_string(upper_limit);

            //std::cout<<region<<std::endl;

	        //generates iterator over region
            //this one is correct:
            hts_itr_t *itr = sam_itr_querys(idx, h, region.c_str());
            /* or do i use this?:
            int bam_fetch(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func);
            */

           //TODO::we need to check for overlapping reads on multiple windows. if this is the case, we need to dynamically increase window size in that case.
            while (sam_itr_next(fp, itr, b) >= 0) {
                per_read_struct* read_output = new per_read_struct();

                //manip reads in region and send to the appropriate functions to get the results
                if (b->core.tid < 0) continue;

                //check if read is in bloom filter, if it is, skip it and go to the next. else add it to the filter
                if(bf_obj.search(bam_get_qname(b)) == 1) 
                    continue;
                else
                    bf_obj.insert(bam_get_qname(b));

                uint32_t read_pos_counter = 0;
                uint32_t ref_pos_counter = 0;
                auto query_sequence_encoded = bam_get_seq(b);
                std::string query_sequence = "";
                int map_quality = int(b->core.qual);

                if(map_quality < opt::min_map_quality) {
                    continue;
                }

                for(int i=0;i<b->core.l_qseq;i++){
                    query_sequence += seq_nt16_str[bam_seqi(bam_get_seq(b), i)];
                }

                //std::cout<<query_sequence<<std::endl;

                auto ref_start_pos = b->core.pos;

                //std::cout<<ref_start_pos<<std::endl;

                /*printf("%s\t%d\t%d\n", h->target_name[b->core.tid], b->core.pos,
                b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)));*/
                //bam handling derived from: https://www.biostars.org/p/4211/
                int x, j, k;
                uint32_t *cigar = bam_get_cigar(b);
                for (k = 0; k < b->core.n_cigar; ++k) {
                    //int op = cigar[k]&16;       //old style htslib
                    //int l = cigar[k]>>4;        //old style htslib
                    int op = bam_cigar_op(cigar[k]);
                    uint32_t l = bam_cigar_oplen(cigar[k]);
                    switch(op) {
                        case BAM_CMATCH:
                            //printf("M");
                            read_pos_counter += l;
                            ref_pos_counter += l;
                            break;
                        case BAM_CHARD_CLIP:
                            //printf("H");
                            break;
                        /*case BAM_CSOFT_CLIP:
                            //printf("S");
                            if(l<opt::min_ins_size){
                                read_pos_counter += l;
                            }
                            else {
                                //call functions here
                                std::cout<<query_sequence<<std::endl;
                                std::cout<<b->core.l_qseq<<" "<<read_pos_counter<<" "<<l<<std::endl;
                                std::string sequence_of_interest = query_sequence.substr(read_pos_counter,l);

                                std::cout<<"softclip: "<<sequence_of_interest<<std::endl;

                                sizing_struct sizing_result;
                                decomposer_struct decomposer_result;
                                methylation_stats methylation_results;

                                decomposer_result = decompose_string(sequence_of_interest,opt::min_repeat_size,opt::max_repeat_size);

                                sizing_result = detect_size(sequence_of_interest,decomposer_result.potential_sequence_in_window);

                                methylation_results = detect_methylation(read_pos_counter,l,b); 

                                if(sizing_result.count<=opt::discovery_sensitivity) {
                                    read_pos_counter += l;
                                    break;
                                }
                                if(decomposer_result.potential_sequence_in_window.empty()) {
                                    read_pos_counter += l;
                                    break;
                                }

                                read_output->region_ref_start = ref_start_pos + ref_pos_counter;
                                read_output->region_ref_end = read_output->region_ref_start + 1;
                                read_output->region_start = read_pos_counter;
                                read_output->region_end = read_pos_counter+l;
                                read_output->interruption_motif = sizing_result.interruption_motif;
                                read_output->size = sizing_result.count;
                                read_output->query_name = bam_get_qname(b);
                                if(bam_is_rev(b)) {
                                    read_output->motif = dna_reverse_complement(decomposer_result.potential_sequence_in_window);
                                    read_output->strand = '-';
                                } 
                                else {
                                    read_output->motif = decomposer_result.potential_sequence_in_window;
                                    read_output->strand = '+';
                                }
                                read_output->min_methylation = methylation_results.min_methylation;
                                read_output->max_methylation = methylation_results.max_methylation;
                                read_output->avg_methylation = methylation_results.avg_methylation;

                                read_pos_counter += l;
                            }
                            break;*/
                        case BAM_CDEL:
                            //printf("D"); does this span positions on the read???
                            ref_pos_counter += l;
                            break;
                        case BAM_CPAD:
                            //printf("P");
                            read_pos_counter += l;
                            break;
                        case BAM_CSOFT_CLIP:
                        case BAM_CREF_SKIP:
                        case BAM_CINS:
                            //printf("I");
                            if(l<opt::min_ins_size){
                                read_pos_counter += l;
                            }
                            else {
                                //call functions here
                                //std::cout<<query_sequence<<std::endl;
                                //std::cout<<b->core.l_qseq<<" "<<read_pos_counter<<" "<<l<<std::endl;
                                std::string sequence_of_interest = query_sequence.substr(read_pos_counter,l);

                                //std::cout<<"insert: "<<sequence_of_interest<<std::endl;

                                sizing_struct sizing_result;
                                decomposer_struct decomposer_result;
                                std::vector<methylation_stats> methylation_results;
                                int haplotype_of_read = 0;

                                decomposer_result = decompose_string(sequence_of_interest,opt::min_repeat_size,opt::max_repeat_size);
                                //TODO::have to filter out homopolymer runs here.
                                if(decomposer_result.potential_sequence_in_window.find_first_not_of(decomposer_result.potential_sequence_in_window[0]) == std::string::npos) {
                                    read_pos_counter += l;
                                    break;
                                }

                                sizing_result = detect_size(sequence_of_interest,decomposer_result.potential_sequence_in_window);

                                /*char *mm_str = bam_aux2Z(bam_aux_get(b, "MM"));
                                char *probability_array = bam_aux2Z(bam_aux_get(b, "ML"));

                                std::cout<<*mm_str<<std::endl;
                                std::cout<<*probability_array<<std::emdl;*/

                                methylation_results = detect_methylation(read_pos_counter,l,b); 

                                if(opt::is_phased) {
                                    haplotype_of_read = get_haplotag(b);
                                }

                                //std::cout<<bam_get_qname(b)<<" "<<decomposer_result.potential_sequence_in_window<<" "<<decomposer_result.potential_count_in_window<<" "<<sizing_result.count<<" "<<sizing_result.interruption_motif<<std::endl;

                                if(sizing_result.count<=opt::discovery_sensitivity) {
                                    read_pos_counter += l;
                                    break;
                                }
                                if(decomposer_result.potential_sequence_in_window.empty()) {
                                    read_pos_counter += l;
                                    break;
                                }

                                read_output->region_ref_start = ref_start_pos + ref_pos_counter;
                                read_output->region_ref_end = read_output->region_ref_start + 1;
                                read_output->region_start = read_pos_counter;
                                read_output->region_end = read_pos_counter+l;
                                read_output->interruption_motif = sizing_result.interruption_motif;
                                read_output->size = sizing_result.count;
                                read_output->query_name = bam_get_qname(b);
                                if(bam_is_rev(b)) {
                                    read_output->motif = dna_reverse_complement(decomposer_result.potential_sequence_in_window);
                                    read_output->strand = '-';
                                } 
                                else {
                                    read_output->motif = decomposer_result.potential_sequence_in_window;
                                    read_output->strand = '+';
                                }
                                read_output->min_methylation = methylation_results[1].min_methylation;
                                read_output->max_methylation = methylation_results[1].max_methylation;
                                read_output->avg_methylation = methylation_results[1].avg_methylation;
                                read_output->up_min_methylation = methylation_results[0].min_methylation;
                                read_output->up_max_methylation = methylation_results[0].max_methylation;
                                read_output->up_avg_methylation = methylation_results[0].avg_methylation;
                                read_output->down_min_methylation = methylation_results[2].min_methylation;
                                read_output->down_max_methylation = methylation_results[2].max_methylation;
                                read_output->down_avg_methylation = methylation_results[2].avg_methylation;
                                read_output->haplotype = haplotype_of_read;

                                read_pos_counter += l;
                            }
                            break;
                        /*case BAM_CREF_SKIP:
                            //printf("S");
                            if(l<opt::min_ins_size){
                                read_pos_counter += l;
                            }
                            else {
                                //call functions here
                                std::cout<<query_sequence<<std::endl;
                                std::cout<<b->core.l_qseq<<" "<<read_pos_counter<<" "<<l<<std::endl;
                                std::string sequence_of_interest = query_sequence.substr(read_pos_counter,l);

                                std::cout<<"refskip: "<<sequence_of_interest<<std::endl;
                                break;

                                sizing_struct sizing_result;
                                decomposer_struct decomposer_result;
                                methylation_stats methylation_results;

                                decomposer_result = decompose_string(sequence_of_interest,opt::min_repeat_size,opt::max_repeat_size);

                                sizing_result = detect_size(sequence_of_interest,decomposer_result.potential_sequence_in_window);

                                methylation_results = detect_methylation(read_pos_counter,l,b); 

                                if(sizing_result.count<=opt::discovery_sensitivity) {
                                    read_pos_counter += l;
                                    break;
                                }
                                if(decomposer_result.potential_sequence_in_window.empty()) {
                                    read_pos_counter += l;
                                    break;
                                }

                                read_output->region_ref_start = ref_start_pos + ref_pos_counter;
                                read_output->region_ref_end = read_output->region_ref_start + 1;
                                read_output->region_start = read_pos_counter;
                                read_output->region_end = read_pos_counter+l;
                                read_output->interruption_motif = sizing_result.interruption_motif;
                                read_output->size = sizing_result.count;
                                read_output->query_name = bam_get_qname(b);
                                if(bam_is_rev(b)) {
                                    read_output->motif = dna_reverse_complement(decomposer_result.potential_sequence_in_window);
                                    read_output->strand = '-';
                                } 
                                else {
                                    read_output->motif = decomposer_result.potential_sequence_in_window;
                                    read_output->strand = '+';
                                }
                                read_output->min_methylation = methylation_results.min_methylation;
                                read_output->max_methylation = methylation_results.max_methylation;
                                read_output->avg_methylation = methylation_results.avg_methylation;

                                read_pos_counter += l;
                            }
                            break;*/
                        default:
                            //printf("?");
                            std::cerr<<"STRr10ToolKit::Cigar_Parse: cannot parse cigar element\n";
                    }
                }
                if(!read_output->motif.empty()) {
                    //TODO::while checking the motifs to add them to the read support also count analogous motifs by getting all possible combinations by rotating the motif. and also their complements. aggregate under only 1 motif
                    //Feature added
                    std::string rotated_motif = read_output->motif;
                    int found = 0;

                    //loop to rotate through the entire string one character at a time and check if they individually exist.
                    for(int i=0; i<read_output->motif.length(); i++) {
                        std::rotate(rotated_motif.begin(), rotated_motif.begin() + 1, rotated_motif.end());
                        if(window_output->motif_aggregate.find(rotated_motif) != window_output->motif_aggregate.end()) {
                            //window_output->motif_aggregate.insert(make_pair(read_output->motif,1));
                            window_output->motif_aggregate[rotated_motif] += 1;
                            found++;
                        }
                        else if(window_output->motif_aggregate.find(dna_reverse_complement(rotated_motif)) != window_output->motif_aggregate.end()) {
                            //window_output->motif_aggregate.insert(make_pair(dna_reverse_complement(read_output->motif),1));
                            window_output->motif_aggregate[dna_reverse_complement(rotated_motif)] += 1;
                            found++;
                            //should i reverse the motif in the read output as well? to merge all into one?
                        }
                    }

                    if(found == 0)
                        window_output->motif_aggregate[read_output->motif] += 1;
                    else
                        found = 0;
                    
                    //--------ORIGINAL CODE START--------------------
                    /*if(window_output->motif_aggregate.find(read_output->motif) != window_output->motif_aggregate.end()) {
                        //window_output->motif_aggregate.insert(make_pair(read_output->motif,1));
                        window_output->motif_aggregate[read_output->motif] += 1;
                    }
                    else if(window_output->motif_aggregate.find(dna_reverse_complement(read_output->motif)) != window_output->motif_aggregate.end()) {
                        //window_output->motif_aggregate.insert(make_pair(dna_reverse_complement(read_output->motif),1));
                        window_output->motif_aggregate[dna_reverse_complement(read_output->motif)] += 1;
                        //should i reverse the motif in the read output as well? to merge all into one?
                    }
                    else {
                        window_output->motif_aggregate[read_output->motif] += 1;
                    }*/
                    //--------ORIGINAL CODE START-------------------
                    window_output->window_aggregate.push_back(*read_output);
                }

                delete read_output;
            }
            int max_read_support=0;
            std::string max_motif;
            for(auto &entry: window_output->motif_aggregate) {
                if(entry.second > max_read_support) {
                    max_read_support = entry.second;
                    max_motif = entry.first;
                }
                //std::cout<<entry.first<<std::endl;
            }
            //TODO::We now have the max motif in the window. we need to aggregate all reads which have that motif and output the aggregate results.
            //feature added
            //FUTURE_FEATURE::check other motifs in motif_aggregate to check if they are similar to the max or if they are analogous to it.
            //feature added above
            if(max_read_support >= opt::min_read_support) {
                //TODO::now that motif has passed the min_read_support test, send all motifs to abPOA to get the consensus and print that out as the motif in the last loop here
                //print out the stats from this window
                //std::vector<per_read_struct> matching_reads; //needed for original code
                int num_of_reads = window_output->window_aggregate.size();
                //char all_motifs[num_of_reads][20];
                uint32_t mean_ref_start = 0;
                uint32_t mean_ref_end = 0;
                int read_count = 0;
		std::vector<std::string> all_motifs;
                std::vector<std::string> consensus_sequences;

                //question here would be, should i send the entire window aggregate to abPOA or should i send only the matching reads?
                //sending the entire window aggregate for now
                for(auto &individual_read: window_output->window_aggregate) {
                    mean_ref_start += individual_read.region_ref_start;
                    mean_ref_end += individual_read.region_ref_end;

		    all_motifs.push_back(individual_read.motif);
                    //strcpy( all_motifs[read_count], individual_read.motif.c_str());

                    read_count += 1;
                }

                mean_ref_end = (uint32_t)(mean_ref_end/read_count);
                mean_ref_start = (uint32_t)(mean_ref_start/read_count);

                consensus_sequences = get_consensus_sequence(all_motifs);

                for(auto &individual_read: window_output->window_aggregate) {
                    std::cout<<individual_read.query_name<<"\t"<<chr<<"\t"<<mean_ref_start<<"\t"<<mean_ref_end<<"\t"<<individual_read.region_start<<"\t"<<individual_read.region_end<<"\t"<<consensus_sequences[0]<<"\t"<<individual_read.interruption_motif<<"\t"<<individual_read.size<<"\t"<<individual_read.haplotype<<"\t"<<individual_read.avg_methylation<<"\t"<<individual_read.min_methylation<<"\t"<<individual_read.max_methylation<<"\t"<<individual_read.up_avg_methylation<<"\t"<<individual_read.up_min_methylation<<"\t"<<individual_read.up_max_methylation<<"\t"<<individual_read.down_avg_methylation<<"\t"<<individual_read.down_min_methylation<<"\t"<<individual_read.down_max_methylation<<std::endl;
                    //chromosome start end reference_length h1_str_length h2_str_length h1_upstream_methylation h1_in_repeat_methylation h1_downstream_methylation h2_upstream_methylation h2_in_repeat_methylation h2_downstream_methylation
                }

                //---------ORIGINAL CODE START---------------------
                /*for(auto &individual_read: window_output->window_aggregate) {
                    if((individual_read.motif == max_motif) || (dna_reverse_complement(individual_read.motif) == max_motif)) {
                        mean_ref_start += individual_read.region_ref_start;
                        mean_ref_end += individual_read.region_ref_end;
                        read_count += 1;

                        matching_reads.push_back(individual_read);
                    }
                }

                mean_ref_end = mean_ref_end/read_count;
                mean_ref_start = mean_ref_start/read_count;

                for(auto matching_read: matching_reads) {
                    std::cout<<matching_read.query_name<<"\t"<<chr<<"\t"<<mean_ref_start<<"\t"<<mean_ref_end<<"\t"<<matching_read.region_start<<"\t"<<matching_read.region_end<<"\t"<<max_motif<<"\t"<<matching_read.interruption_motif<<"\t"<<matching_read.size<<"\t"<<matching_read.avg_methylation<<"\t"<<matching_read.min_methylation<<"\t"<<matching_read.max_methylation<<"\t"<<matching_read.haplotype<<std::endl;
                    //chromosome start end reference_length h1_str_length h2_str_length h1_upstream_methylation h1_in_repeat_methylation h1_downstream_methylation h2_upstream_methylation h2_in_repeat_methylation h2_downstream_methylation
                }*/
                //--------ORIGINAL CODE END-------------------------------
            }
            //otherwise, there aren't any STRs that pass all the filters. moving to the next window

            delete window_output;
	        lower_limit = upper_limit;
	        upper_limit += opt::window_size;
            hts_itr_destroy(itr);
        }
    }

    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fp);
}
