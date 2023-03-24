#include "./str_utilities.h"

// for nt
// AaCcGgTtNn ==> 0,1,2,3,4
unsigned char nt4_table[256] = {
       0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
const char nt256_table[256] = {
       'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
};

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

uint32_t get_median(std::vector<uint32_t> coordinates_array) {
    sort(coordinates_array.begin(), coordinates_array.end());
    int size = coordinates_array.size();

    if (coordinates_array.size() % 2 != 0)
        return (uint32_t)coordinates_array[size/2];
    else
        return (uint32_t)(coordinates_array[(size-1)/2] + coordinates_array[size/2])/2.0;
}

sizing_struct detect_size(std::string sequence_of_interest, std::string potential_str_sequence) {
    int rindex = 0;
    int lindex = 0;

    int motif_length = potential_str_sequence.length();
    int count = 0;
    std::string int_motif = "";
    int int_count = 0;

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
            if(int_count == 0) {
                if((rindex-(lindex+motif_length)) <= 50) {
                    int_motif = sequence_of_interest.substr(lindex+motif_length,rindex-(lindex+motif_length)); //the addition in first parameter is to get past the existing motif at the left index, to the start of the interruption. 
                    count++;
                }
                int_count++;
            }
            int_count++;
            if(int_count > 2)
                break;
        }

        lindex = rindex;
        ++rindex;
    }

    sizing_struct return_variable = {count+1,int_motif,int_count}; //the +1 compensates for the initial skipped find

    return return_variable;
}

//Legacy method to manually parse Methylation tags
/*methylation_stats detect_methylation(int region_start, int region_end, bam1_t *b) { //pass bam record here    
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
}*/

//legacy decomposer. trying to refactor this in the new function
/*decomposer_struct decompose_string(std::string sequence_of_interest, int lower_limit, int upper_limit) {

    std::map<std::string, int> subsequences;
    sizing_struct sizing_result;

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
}*/

decomposer_struct decompose_string(std::string sequence_of_interest, int lower_limit, int upper_limit) {

    std::map<std::string, int> subsequences;
    sizing_struct sizing_result;

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
    int insert_span = 0;
    int overall_max_value=0;
    std::string overall_max_key;

    for(int motif_length=lower_limit ; motif_length <= upper_limit ; motif_length++) {
        for(auto &entry: subsequences) {
            if(entry.second > max_value && entry.first == motif_length) {
                max_key = entry.first;
                max_value = entry.second;
            }
        }
        int tmp_span = (int)((max_value * max_key.length() * 100) / sequence_of_interest.length());
        if(tmp_span > insert_span) {
            insert_span = tmp_span;
            overall_max_value = max_value;
            overall_max_key = max_key;
        }
    }

    decomposer_struct return_variable = {overall_max_key,overall_max_value};

    return return_variable;
}

//new decomposer algorithm
/*decomposer_struct decompose_string(std::string sequence_of_interest, int lower_limit, int upper_limit) {

    std::vector<std::map<int, std::string>> all_k_mers;
    std::vector<std::tuple<std::string,int> k_mer_max_detected;

    //getting all the k-mers for the different values of lower_limit<=k<=upper_limit
    for(int motif_length=lower_limit ; motif_length <= upper_limit ; motif_length++) {
        int lower_window_var = 0;
        int upper_window_var = motif_length;

        std::map<int,std::string> k_mer_index;

        while(upper_window_var <= sequence_of_interest.length()) {
            std::string k_mer = sequence_of_interest.substr(lower_window_var,motif_length);

            k_mer_index.insert({lower_window_var,k_mer});
            
            lower_window_var++;
            upper_window_var++;
        }

        all_k_mers.push_back(k_mer_index);
    }

    int max_value=0;
    std::string max_key;

    for(auto &entry: subsequences) {
        if(entry.second > max_value) {
            max_key = entry.first;
            max_value = entry.second;
        }
    }

    for(auto kmer_index_map: all_k_mers) {
        //go through index by index to get max per motif length
    }

    decomposer_struct return_variable = {max_key,max_value};

    return return_variable;
}*/

int get_haplotag(bam1_t *b) {
    int haplotag = bam_aux2i(bam_aux_get(b, "HP")); //encoded as HP:i:1 or HP:i:2 by Whatshap

    if(haplotag == 1 || haplotag == 2) {
        return haplotag;
    }
    else {
        return 0;
    }
}

std::vector<std::string> get_consensus_sequence(std::vector<std::string> sequences) { //m is rows and n is columns
    int i, j, n_seqs = sequences.size();

    std::vector<std::string> output_msa;

    // initialize variables
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

    // alignment parameters
    // abpt->align_mode = 0; // 0:global 1:local, 2:extension
    // abpt->mat_fn = strdup("HOXD70.mtx"); abpt->use_score_matrix = 1; // score matrix instead of constant match/mismatch score
    // abpt->match = 2;      // match score
    // abpt->mismatch = 4;   // mismatch penalty
    // abpt->gap_mode = ABPOA_CONVEX_GAP; // gap penalty mode
    // abpt->gap_open1 = 4;  // gap open penalty #1
    // abpt->gap_ext1 = 2;   // gap extension penalty #1
    // abpt->gap_open2 = 24; // gap open penalty #2
    // abpt->gap_ext2 = 1;   // gap extension penalty #2
                             // gap_penalty = min{gap_open1 + gap_len * gap_ext1, gap_open2 + gap_len * gap_ext2}
    // abpt->bw = 10;        // extra band used in adaptive banded DP
    // abpt->bf = 0.01; 
     
    // output options
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
    abpt->w = 6, abpt->k = 9; abpt->min_w = 10; // minimizer-based seeding and partition
    //abpt->progressive_poa = 1;
    abpt->max_n_cons = 1; // to generate 2 consensus sequences

    abpoa_post_set_para(abpt);

    // collect sequence length, trasform ACGT to 0123
    int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
    //int **weights = (int**)malloc(sizeof(int*) * n_seqs);
    for (i = 0; i < n_seqs; ++i) {
        //seq_lens[i] = strlen(sequences[i]);
	seq_lens[i] = sequences[i].length();
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
        //weights[i] = (int*)malloc(sizeof(int) * seq_lens[i]);
        for (j = 0; j < seq_lens[i]; ++j) {
            bseqs[i][j] = nt4_table[(int)sequences[i][j]];
            //if (j >= 12) weights[i][j] = 2;
            //else weights[i][j] = 0;
        }
    }

    // 1. directly output to stdout
    //fprintf(stdout, "=== output to stdout ===\n");
    //abpt->use_qv = 1;
    // perform abpoa-msa
    // set weights as NULL if no quality score weights are used
    //abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, weights, stdout);
    //std::cout<<"before msa\n";
    int val = abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL);
    //std::cout<<"after msa\n";
    if(val == 1) {
        //free seq-related variables
	for (i = 0; i < n_seqs; ++i) { 
		free(bseqs[i]); 
		//free(weights[i]); 
	}
	free(bseqs); 
	free(seq_lens); 
	//free(weights);
	
	// free abpoa-related variables
	abpoa_free(ab); abpoa_free_para(abpt);
	output_msa.push_back("");
    }

    // 2. output MSA alignment and consensus sequence stored in (abpoa_cons_t *)
    abpoa_cons_t *abc = ab->abc;
    //fprintf(stdout, "=== stored in variables ===\n");
    //fprintf(stdout, ">Multiple_sequence_alignment\n");
    /*for (i = 0; i < abc->n_seq; ++i) {
        for (j = 0; j < abc->msa_len; ++j) {
            fprintf(stdout, "%c", nt256_table[abc->msa_base[i][j]]);
        }
        fprintf(stdout, "\n");
    }*/

    for (i = 0; i < abc->n_cons; ++i) {
        std::string cons_seq;
        //fprintf(stdout, ">Consensus_sequence");
        /*if (abc->n_cons > 1) {
            fprintf(stdout, "_%d ", i+1);
            for (j = 0; j < abc->clu_n_seq[i]; ++j) { // output read ids for each cluster/group
                fprintf(stdout, "%d", abc->clu_read_ids[i][j]);
                if (j != abc->clu_n_seq[i]-1) fprintf(stdout, ",");
            }
        }*/
        //fprintf(stdout, "\n");
        for (j = 0; j < abc->cons_len[i]; ++j) {
            //fprintf(stdout, "%c", nt256_table[abc->cons_base[i][j]]);
            cons_seq += nt256_table[abc->cons_base[i][j]];
        }
        //fprintf(stdout, "\n");
        output_msa.push_back(cons_seq);
    }

    /* generate DOT partial order graph plot */
    //abpt->out_pog = strdup("example.png"); // dump parital order graph to file
    //if (abpt->out_pog != NULL) abpoa_dump_pog(ab, abpt);

    // free seq-related variables
    for (i = 0; i < n_seqs; ++i) { 
	    free(bseqs[i]); 
	    //free(weights[i]); 
    }
    free(bseqs); 
    free(seq_lens); 
    //free(weights);

    // free abpoa-related variables
    abpoa_free(ab); abpoa_free_para(abpt); 

    return output_msa;
}

std::vector<methylation_stats> detect_methylation(int region_start, int region_end, bam1_t *b) {
    float min_methylation[3] = {0,0,0};
    float max_methylation[3] = {0,0,0};
    float avg_methylation[3] = {0,0,0};
    float total_methylation[3] = {0,0,0};

    std::vector<methylation_stats> return_variable;

    int read_pos_count = 0;

    hts_base_mod_state* base_mod_states = hts_base_mod_state_alloc();
    const int max_mods = 5;
    hts_base_mod mods[max_mods];

    int initialize_base_mod_states = bam_parse_basemod(b, base_mod_states);

    if(initialize_base_mod_states == -1) {
        for(int k = 0; k < 3; k++)
            return_variable.push_back({0,0,0});

        return return_variable;
    }

    int out_position;

    while(bam_next_basemod(b, base_mod_states, mods, max_mods, &out_position) > 0) { //iterating over the number of base mods found
        read_pos_count += out_position + 1;

        for(int i = 0; i < max_mods; i++) {
            if(mods[i].modified_base == 'm') {

                float probability_of_mod = mods[i].qual != -1 ? ((float)mods[i].qual/(float)256) : 0;

                if (read_pos_count >= region_start && read_pos_count <= region_end) { //check if basemod is in the insert region
                    if (probability_of_mod > max_methylation[1]) {
                        max_methylation[1] = probability_of_mod;
                    }
                    if (probability_of_mod < min_methylation[1]) {
                        min_methylation[1] = probability_of_mod;
                    }
                    total_methylation[1] += probability_of_mod;
                }
                if (read_pos_count < region_start && read_pos_count >= (region_start-3000)) { //check if basemod is in the upstream region
                    if (probability_of_mod > max_methylation[0]) {
                        max_methylation[0] = probability_of_mod;
                    }
                    if (probability_of_mod < min_methylation[0]) {
                        min_methylation[0] = probability_of_mod;
                    }
                    total_methylation[0] += probability_of_mod;
                }
                if (read_pos_count > region_end && read_pos_count <= (region_end+3000)) { //check if basemod is in the downstream region
                    if (probability_of_mod > max_methylation[2]) {
                        max_methylation[2] = probability_of_mod;
                    }
                    if (probability_of_mod < min_methylation[2]) {
                        min_methylation[2] = probability_of_mod;
                    }
                    total_methylation[2] += probability_of_mod;
                }
            }
        }
        
    }
    //in-region avg_methylation calculation
    if (read_pos_count >= region_end)
        avg_methylation[1] = total_methylation[1] / (float)(region_end - region_start);
    else if (read_pos_count < region_end && read_pos_count >= region_start)
        avg_methylation[1] = total_methylation[1] / (float)(read_pos_count - region_start);
    else if (read_pos_count < region_start)
        avg_methylation[1] = 0;

    //downstream avg_methylation calculation
    if(read_pos_count >= (region_end+3000))
        avg_methylation[2] = total_methylation[2] / (float)(3000);
    else if(read_pos_count < (region_end+3000) && read_pos_count >= region_end)
        avg_methylation[2] = total_methylation[2] / (float)(read_pos_count - region_end);
    else if(read_pos_count < region_end)
        avg_methylation[2] = 0;

    //upstream avg_methylation calculation
    if(read_pos_count >= (region_start))
        avg_methylation[0] = total_methylation[0] / (float)((region_start - 3000 >= 0) ? 3000 : region_start);
    else if(read_pos_count < region_start && read_pos_count >= (region_start - 3000))
        avg_methylation[0] = total_methylation[0] / (float)((region_start - 3000 >= 0) ? (read_pos_count - (region_start - 3000)): read_pos_count);
    else if(read_pos_count < (region_start - 3000))
        avg_methylation[0] = 0;

    //methylation_stats return_variable = {max_methylation,min_methylation,avg_methylation};

    for(int j = 0; j < 3; j++)
        return_variable.push_back({max_methylation[j],min_methylation[j],avg_methylation[j]});

    hts_base_mod_state_free(base_mod_states);

    //delete methylation_prob;

    return return_variable;
}
