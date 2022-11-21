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

struct decomposer_struct {
    std::string potential_sequences_in_window;
    int potential_count_in_window;
}decomposer_output;

struct sizing_struct {
    int count;
    std::string interruption_motif;
}sizing_out;

struct methylation_stats {
    float max_methylation;
    float min_methylation;
    float avg_methylation;
};

/*std::vector<size_t> get_all_found_pos(std::string to_search_in, std::string to_search_for) {
    std::vector<size_t> positions; // holds all the positions that sub occurs within str

    size_t pos = str.find(sub, 0);
    while(pos != string::npos)
    {
        positions.push_back(pos);
        pos = str.find(sub,pos+1);
    }

    return positions;
}*/

std::tuple<int, std::string> detect_size(std::string sequence_of_interest, std::string potential_str_sequence) {
    int rindex = 0;
    int lindex = 0;

    int motif_length = potential_str_sequence.length();
    sizing_struct* output_variable = new sizing_struct();

    while ((rindex = sequence_of_interest.find(potential_str_sequence, rindex)) != std::string::npos) {
        std::cout << rindex << std::endl;
        if(rindex - lindex == motif_length) {
            output_variable->count++;
        }
        else {
            output_variable->interruption_motif = potential_str_sequence.substr(lindex,rindex);
        }

        lindex = rindex;
        rindex += motif_length;
    }

    return {output_variable->count,output_variable->interruption_motif};
}

std::tuple<float, float, float> detect_methylation(int region_start, int region_end, bam1_t *b) { //pass bam record here
    
    float min_methylation = 0;
    float max_methylation = 0;
    float avg_methylation = 0;
    float total_methylation = 0;

    int read_pos_count = 0;
    
    char *mm_str = bam_aux2Z(bam_aux_get(b, "MM"));
    char *probability_array = bam_aux2Z(bam_aux_get(b, "ML"));

    std::string tmp_mm1(mm_str);
    std::string tmp_prob(probability_array);

    //Multiple mods are not supported as of now
    if(std::count(tmp_mm1.begin(),tmp_prob.end(),';') > 1) {
        return {0,0,0};
    }

    std::vector<std::string> positions;
    std::vector<std::string> probabilites;

    std::string token;

    std::stringstream tmp_mm_str(mm_str);
    while (getline(tmp_mm_str, token, ',')){
        positions.push_back(token);
    }

    std::stringstream tmp_probability_array(probability_array);
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

    return {min_methylation,max_methylation,avg_methylation};
}

/*std::tuple<*methylation_stats, *methylation_stats, *methylation_stats> detect_methylation_all_regions(int read_start, int read_end, bam1_t *b) { //pass bam record here
    
    float min_methylation = 0;
    float max_methylation = 0;
    float avg_methylation = 0;
    float total_methylation = 0;

    int read_pos_count = 0;
    
    char *mm_str = bam_aux2Z(bam_aux_get_core(b, "MM"));
    char *probability_array = bam_aux2Z(bam_aux_get_core(b, "ML"));

    //Multiple mods are not supported as of now
    if(std::count(mm_str.begin(),mm_str.end(),';') > 1) {
        return {0,0,0};
    }

    std::vector<std::string> positions;
    std::vector<std::string> probabilites;

    std::string token;

    std::stringstream tmp_mm_str(mm_str);
    while (getline(tmp_mm_str, token, ',')){
        positions.push_back(token);
    }

    std::stringstream tmp_probability_array(probability_array);
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

        if read_pos_count >= read_start && read_pos_count <= read_end {
            if probability_of_mod > max_methylation {
                max_methylation = methylation_base_probability;
            }
            if probability_of_mod < min_methylation {
                min_methylation = methylation_base_probability;
            }
            total_methylation += probability_of_mod;
        }

        ++mm_iterator;
        ++ml_iterator;
    }

    if read_pos_count > read_end {
        avg_methylation = total_methylation / (float)(read_end - read_start);
    }
    else if read_pos_count < read_end && read_pos_count > read_start {
        avg_methylation = total_methylation / (float)(read_pos_count - read_start);
    }
    else if read_pos_count < read_start {
        avg_methylation = 0;
    }

    return {min_methylation,max_methylation,avg_methylation};
}*/

std::tuple<std::string, int> decompose_string(std::string sequence_of_interest, int lower_limit, int upper_limit) {

    std::map<std::string, int> subsequences;

    for(int motif_length=lower_limit ; motif_length <= upper_limit ; motif_length++) {
        int lower_window_var = 0;
        int upper_window_var = motif_length;

        while(upper_window_var <= sequence_of_interest.length()) {
            std::string sequence_in_window = sequence_of_interest.substr(lower_window_var,upper_window_var);

            if (subsequences.find(sequence_in_window) == subsequences.end()) {
                subsequences.insert({sequence_in_window, 1});
            }
            else {
                subsequences[sequence_in_window]+=1;
            }
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

    return {max_key,max_value};
}