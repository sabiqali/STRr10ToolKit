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

std::vector<std::string> get_consensus_sequence(int m, int n, char sequences[][20]) { //m is rows and n is columns
    int i, j, n_seqs = m;

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
    abpt->progressive_poa = 1;
    abpt->max_n_cons = 2; // to generate 2 consensus sequences

    abpoa_post_set_para(abpt);

    // collect sequence length, trasform ACGT to 0123
    int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
    int **weights = (int**)malloc(sizeof(int*) * n_seqs);
    for (i = 0; i < n_seqs; ++i) {
        seq_lens[i] = strlen(sequences[i]);
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
        weights[i] = (int*)malloc(sizeof(int) * seq_lens[i]);
        for (j = 0; j < seq_lens[i]; ++j) {
            bseqs[i][j] = nt4_table[(int)sequences[i][j]];
            if (j >= 12) weights[i][j] = 2;
            else weights[i][j] = 0;
        }
    }

    // 1. directly output to stdout
    fprintf(stdout, "=== output to stdout ===\n");
    abpt->use_qv = 1;
    // perform abpoa-msa
    // set weights as NULL if no quality score weights are used
    //abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, weights, stdout);
    int val = abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL);

    // 2. output MSA alignment and consensus sequence stored in (abpoa_cons_t *)
    abpoa_cons_t *abc = ab->abc;
    fprintf(stdout, "=== stored in variables ===\n");
    fprintf(stdout, ">Multiple_sequence_alignment\n");
    for (i = 0; i < abc->n_seq; ++i) {
        for (j = 0; j < abc->msa_len; ++j) {
            fprintf(stdout, "%c", nt256_table[abc->msa_base[i][j]]);
        }
        fprintf(stdout, "\n");
    }

    for (i = 0; i < abc->n_cons; ++i) {
        std::string cons_seq;
        fprintf(stdout, ">Consensus_sequence");
        if (abc->n_cons > 1) {
            fprintf(stdout, "_%d ", i+1);
            for (j = 0; j < abc->clu_n_seq[i]; ++j) { // output read ids for each cluster/group
                fprintf(stdout, "%d", abc->clu_read_ids[i][j]);
                if (j != abc->clu_n_seq[i]-1) fprintf(stdout, ",");
            }
        }
        fprintf(stdout, "\n");
        for (j = 0; j < abc->cons_len[i]; ++j) {
            fprintf(stdout, "%c", nt256_table[abc->cons_base[i][j]]);
            cons_seq += nt256_table[abc->cons_base[i][j]];
        }
        fprintf(stdout, "\n");
        output_msa.push_back(cons_seq);
    }

    /* generate DOT partial order graph plot */
    abpt->out_pog = strdup("example.png"); // dump parital order graph to file
    if (abpt->out_pog != NULL) abpoa_dump_pog(ab, abpt);

    // free seq-related variables
    for (i = 0; i < n_seqs; ++i) { free(bseqs[i]); free(weights[i]); }
    free(bseqs); free(seq_lens); free(weights);

    // free abpoa-related variables
    abpoa_free(ab); abpoa_free_para(abpt); 

    return output_msa;
}