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
#include "htslib/sam.h"
#include "htslib/bam.h"

std::tuple<int, std::string> detect_size(std::string sequence_of_interest, std::string potential_str_sequence);

std::tuple<float, float, float> detect_methylation(int read_start, int read_end, bam1_t *b);

std::tuple<std::string, int> decompose_string(std::string sequence_of_interest, int lower_limit, int upper_limit);

std::string dna_reverse_complement(std::string seq);

#endif

