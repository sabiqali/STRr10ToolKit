#include "./str_utilities.h"

std::string dna_reverse_complement(std::string seq) {
    reverse(seq.begin(),seq.end());
    for (std::size_t i = 0; i < seq.length(); ++i){
        switch (seq[i]){
        case 'A':
            seq[i] = 'T';
            break;
        case 'C':
            seq[i] = 'G';
            break;
        case 'G':
            seq[i] = 'C';
            break;
        case 'T':
            seq[i] = 'A';
            break;
        }
    }
    return seq;
}

sizing_struct detect_size(std::string sequence_of_interest, std::string potential_str_sequence) {
    int rindex = 0;
    int lindex = 0;

    int motif_length = potential_str_sequence.length();
    int count = 0;
    std::string int_motif;

    while ((rindex = sequence_of_interest.find(potential_str_sequence, rindex)) != std::string::npos) { //could implement KMP algorithm here to get indeces in O(m)
        if(lindex == 0) {
            lindex = rindex;
            ++rindex;
            continue;
        }
        if((rindex - lindex) == motif_length) {
            count++;
        }
        else {
            int_motif = sequence_of_interest.substr(lindex+motif_length,rindex-(lindex+motif_length)); //the addition in first parameter is to get past the existing motif at the left index, to the start of the interruption. 
            count++;
        }

        lindex = rindex;
        ++rindex;
    }

    sizing_struct return_variable = {count+1,int_motif}; //the +1 compensates for the initial skipped find

    return return_variable;
}

//methylation_stats detect_methylation(int region_start, int region_end, std::string mm_string, std::string prob_array) { //pass bam record here
methylation_stats detect_methylation(int region_start, int region_end, bam1_t *b) { //pass bam record here    
    float min_methylation = 0;
    float max_methylation = 0;
    float avg_methylation = 0;
    float total_methylation = 0;

    int read_pos_count = 0;
    
    char *mm_str = bam_aux2Z(bam_aux_get(b, "MM"));
    char *probability_array = bam_aux2Z(bam_aux_get(b, "ML"));

    //char *mm_str = const_cast<char*>(mm_string.c_str());
    //char *probability_array = const_cast<char*>(prob_array.c_str());

    std::string tmp_mm1(mm_str);
    std::string tmp_prob(probability_array);

    //Multiple mods are not supported as of now
    if(std::count(tmp_mm1.begin(),tmp_prob.end(),';') > 1) {
        //return std::make_tuple(0,0,0);
        methylation_stats return_variable = {0,0,0};

        return return_variable;
    }

    std::vector<std::string> positions;
    std::vector<std::string> probabilites;

    std::string token;

    tmp_mm1 = tmp_mm1.substr(0, tmp_mm1.find(';'));

    std::stringstream tmp_mm_str(tmp_mm1);
    while (getline(tmp_mm_str, token, ',')){
        positions.push_back(token);
    }

    std::stringstream tmp_probability_array(tmp_prob);
    while (getline(tmp_probability_array, token, ',')){
        probabilites.push_back(token);
    }

    std::vector<std::string>::iterator mm_iterator = positions.begin();
    std::advance(mm_iterator, 1);

    std::vector<std::string>::iterator ml_iterator = probabilites.begin();
    std::advance(ml_iterator, 1);

    while (mm_iterator != positions.end())
    {
        int position_of_mod = std::stoi(*mm_iterator);
        float probability_of_mod = std::stof(*ml_iterator)/(float)255;

        read_pos_count += position_of_mod + 1;

        if (read_pos_count >= region_start && read_pos_count <= region_end) {
            if (probability_of_mod > max_methylation) {
                max_methylation = probability_of_mod;
            }
            if (probability_of_mod < min_methylation) {
                min_methylation = probability_of_mod;
            }
            total_methylation += probability_of_mod;
        }

        ++mm_iterator;
        ++ml_iterator;
    }

    if (read_pos_count > region_end) {
        avg_methylation = total_methylation / (float)(region_end - region_start);
    }
    else if (read_pos_count < region_end && read_pos_count > region_start) {
        avg_methylation = total_methylation / (float)(read_pos_count - region_start);
    }
    else if (read_pos_count < region_start) {
        avg_methylation = 0;
    }

    methylation_stats return_variable = {max_methylation,min_methylation,avg_methylation};

    return return_variable;
}

decomposer_struct decompose_string(std::string sequence_of_interest, int lower_limit, int upper_limit) {

    std::map<std::string, int> subsequences;

    for(int motif_length=lower_limit ; motif_length <= upper_limit ; motif_length++) {
        int lower_window_var = 0;
        int upper_window_var = motif_length;

        while(upper_window_var <= sequence_of_interest.length()) {
            std::string sequence_in_window = sequence_of_interest.substr(lower_window_var,motif_length);

            if (subsequences.find(sequence_in_window) == subsequences.end()) {
                subsequences.insert(make_pair(sequence_in_window, 1));
            }
            else {
                subsequences[sequence_in_window]+=1;
            }

            lower_window_var++;
            upper_window_var++;
        }
    }

    int max_value=0;
    std::string max_key;

    for(auto &entry: subsequences) {
        if(entry.second > max_value) {
            max_key = entry.first;
            max_value = entry.second;
        }
    }

    decomposer_struct return_variable = {max_key,max_value};

    return return_variable;
}

int get_haplotag(bam1_t *b) {
    int haplotag = bam_aux2i(bam_aux_get(b, "HP")); //encoded as HP:i:1 or HP:i:2 by Whatshap

    if(haplotag == 1 || haplotag == 2) {
        return haplotag;
    }
    else {
        return 0;
    }
}
