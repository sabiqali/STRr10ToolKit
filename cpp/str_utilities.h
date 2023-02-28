#ifndef STR_UTILITIES_H
#define STR_UTILITIES_H

#include <zlib.h>
#include <stdio.h>
#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <stdbool.h>
#include <fstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <unordered_set>
#include <chrono>
#include <sstream>
#include <getopt.h>
#include <tuple>
//#include "./htslib/bam.h"
#include "htslib/htslib/vcf.h"
#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"
#include "htslib/htslib/hts.h"
#include <stdint.h>
#include "abPOA/include/abpoa.h"

struct decomposer_struct {
    std::string potential_sequence_in_window;
    int potential_count_in_window;
};

struct sizing_struct {
    int count;
    std::string interruption_motif;
};

struct methylation_stats {
    float max_methylation;
    float min_methylation;
    float avg_methylation;
};

//std::tuple<int, std::string> detect_size(std::string sequence_of_interest, std::string potential_str_sequence);

//std::tuple<float, float, float> detect_methylation(int read_start, int read_end, bam1_t *b);

//std::tuple<std::string, int> decompose_string(std::string sequence_of_interest, int lower_limit, int upper_limit);

sizing_struct detect_size(std::string sequence_of_interest, std::string potential_str_sequence);

//methylation_stats detect_methylation(int read_start, int read_end, bam1_t *b); //method to do the methylation tags retrieval manually

methylation_stats detect_methylation(int region_start, int region_end, bam1_t *b); //method to retrieve methylation tags using the htslib API: https://github.com/samtools/htslib/blob/develop/htslib/sam.h#L2260

decomposer_struct decompose_string(std::string sequence_of_interest, int lower_limit, int upper_limit);

std::string dna_reverse_complement(std::string seq);

int get_haplotag(bam1_t *b);

std::vector<std::string> get_consensus_sequence(std::vector<std::string> sequences);  //m is rows and n is columns

#endif

