//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Sabiq Chaudhary (schaudhary@oicr.on.ca)
//---------------------------------------------------------

#[macro_use]
extern crate rust_htslib;
extern crate substring;
use crate::str_sizing::size_struct;

use rust_htslib::{bam, faidx, bam::Read, bam::record::CigarStringView};
use substring::Substring;
use std::cmp;
use std::collections::HashMap;

mod methylation_detection;
mod str_sizing;

struct per_read_struct {
    pub read_start: u32,
    pub read_end: u32,
    pub motif: &str,
    pub haplotype: u32,
    pub size: u32,
    pub min_methylation: u32,
    pub max_methylation: u32,
    pub avg_methylation: u32,
    pub interruption_motif: &str,
}

struct per_window_struct {
    pub read_name: &str,
    pub window_start: u32,
    pub window_end: u32,
    pub motif: &str,
    pub haplotype: u32,
    pub size: u32,
    pub min_methylation: u32,
    pub max_methylation: u32,
    pub avg_methylation: u32,
    pub interruption_motif: &str,
}

pub struct all_charac_struct {
    pub found_flag: bool,
    pub window_start: u32,
    pub window_end: u32,
    pub motif: &str,
    pub haplotype: u32,
    pub size: u32,
    pub min_methylation: u32,
    pub max_methylation: u32,
    pub avg_methylation: u32,
    pub interruption_motif: u32,
    pub read_support: u32,
}

pub struct decomposer_struct {
    pub potential_sequences_in_window: Vec<&str>,
    pub potential_count_in_window: Vec<u32>,
}

/// Returns the length of the longest substring without repeating characters.
///
/// # Arguments
///
/// * 's' - The input string
pub fn length_of_longest_substring(s: String) -> 132 {
    let mut length = 0;
    
    let mut char_set: HashMap<char, usize> = HashMap::new();
    
    let mut start = 0:
    for (end, c) in s.char_indices() {
        if let Some(&n) = char_set.get(&c) {
            start = cmp::max(start, n);
        }
        length = cmp::max(length, end - start + 1);
        char_set.insert(c, end + 1);
    }

    length as i32
}

pub fn decompose_string(sequence_of_interest: &str, lower_limit: u32, upper_limit: u32) -> (&str,u32) {
    //TODO::take in the string of interest from the CIGAR String and then find if there is a repeating unit present

    let mut decomposer_return_struct: decomposer_struct;

    let mut subsequences: HashMap<&str, u32> = HashMap::new();

    for motif_length in (lower_limit,upper_limit) {

        let mut lower_window_var = 0;
        let mut upper_window_var = motif_length

        while upper_window_var <= sequence_of_interest.len() {
            sequence_in_window = sequence_of_interest.substring(lower_window_var,upper_window_var);

            //let occurance = subsequences.entry(sequence_in_window).or_insert(1);

            if subsequences.contains_key(sequence_in_window) {
                //let number_of_occurances = subsequences.get(sequence_in_window);
                *subsequences.entry(sequence_in_window).or_insert(0) += 1;
            }
            else {
                subsequences.insert(sequence_in_window,1);
            }
        }
    }

    let mut max_key: &str;
    let mut max_val = 0; //can we use this as count instead of the sizing function?? since we expect reads to be near perfect. NOTE::this does not take into account interruptions.

    for (key, value) in &*subsequences {
        if value > max_val {
            max_key = key;
            max_val = value;
        }
    }

    return (max_key,max_val);
}

pub fn detect_loci(window_start: u32, window_end: u32, alignments: rust_htslib::bam::pileup::Alignments<'_>, lower_limit: u32, upper_limit: u32, min_read_support: u32, min_ins_size: u32, discovery_sensitivity: u32) -> per_window_struct {
    println!("Inside STR Discovery");
    //also pass a previous window ref_start and ref_end parameter. if present calculation aligns to the same region, skip this window and move to the next.

    let mut window_aggregate: Vec<per_window_struct> = Vec::new();

    //let mut read_support_counter: u32;

    for a in alignments {
        if a.record().seq().len() == 0 {
            continue;
        }

        //TODO::check bam record for any repeats here
        let mut cigar_string = a.cigar();
        if cigar_string.trim().is_empty() {
            continue;
        }

        let mut read_characteristics: per_read_struct;
        
        let mut event_length = 0; //this is present to count the number of bases in each 'M','I','S', or other CIGAR events
        let mut read_counter = 0; //this is present to keep track of number of aligned bases till we hit the region of interest in the CIGAR string
        for letter in cigar_string {
            if letter.is_numeric() {
                event_length = 10*event_length + letter.to_digit(10).unwrap();
            }
            else {
                match letter {
                    'I' => {
                        if event_length > min_ins_size {
                            let mut sequence_of_interest = a.record().seq().substring(read_counter,event_length); //slice it from (read_counter, event_length)
                            //Do we want to get complex motifs? if so, look for multiple expanded motifs in the window and then try to size both in the window
                            let mut (potential_str_sequence,potential_count_from_discovery) = decompose_string(&sequence_of_interest,lower_limit,upper_limit);
                            if potential_str_sequence.trim().is_empty() {
                                continue;
                            }
                            //TODO::call sizing and methylation function if there is a potential str, to calc characteristics
                            //TODO::send sequence of interest as parameter to the sizing
                            //TODO::send record with read start and end positions to calculate methylation
                            //FORMAT::store in per read map which is potential_sequence maps to {ref_start,ref_end,motif,haplotype,size,methylation,interruption_motif}
                            
                            let mut size_func_return: size_struct;
                            size_func_return = detect_loci(&sequence_of_interest,&potential_str_sequence);

                            if size_func_return.count <= discovery_sensitivity {
                                continue;
                            }

                            (avg_methylation,min_methylation, max_methylation) = detect_methylation(read_counter, read_counter+event_length, a.record())

                            read_characteristics.read_start = read_counter;
                            read_characteristics.read_end = read_counter+event_length;
                            read_characteristics.motif = potential_str_sequence;
                            read_characteristics.size = size_func_return.count;
                            read_characteristics.avg_methylation = avg_methylation;
                            read_characteristics.min_methylation = min_methylation;
                            read_characteristics.max_methylation = max_methylation;

                            read_counter += event_length;
                            event_length = 0;
                        }
                    },
                    'S' => {
                        if event_length > min_ins_size {
                            let mut sequence_of_interest = a.record().seq().substring(read_counter,event_length); //slice it from (read_counter, event_length)
                            let mut (potential_str_sequence,potential_count_from_discovery) = decompose_string(&sequence_of_interest,lower_limit,upper_limit);
                            if potential_str_sequence.trim().is_empty() {
                                continue;
                            }
                            //TODO::call sizing and methylation function if there is a potential str, to calc characteristics
                            //TODO::send sequence of interest as parameter to the sizing
                            //TODO::send record with read start and end positions to calculate methylation
                            //FORMAT::store in per read map which is read_name maps to {potential_motif,ref_start,ref_end,motif,haplotype,size,methylation,interruption_motif}
                            let mut size_func_return: size_struct;
                            size_func_return = detect_loci(&sequence_of_interest,&potential_str_sequence);

                            if size_func_return.count <= discovery_sensitivity {
                                continue;
                            }

                            (avg_methylation,min_methylation, max_methylation) = detect_methylation(read_counter, read_counter+event_length, a.record())
                            
                            read_characteristics.read_start = read_counter;
                            read_characteristics.read_end = read_counter+event_length;
                            read_characteristics.motif = potential_str_sequence;
                            read_characteristics.size = size_func_return.count;
                            read_characteristics.avg_methylation = avg_methylation;
                            read_characteristics.min_methylation = min_methylation;
                            read_characteristics.max_methylation = max_methylation;

                            read_counter += event_length;
                            event_length = 0;
                        }
                    },
                    _ => {
                        read_counter += event_length;
                        event_length = 0;
                        continue;    
                    },
                }
            }
        }

        //TODO::aggregate all under per window map which is read_name maps to {potential_motif,ref_start,ref_end,window_start,window_end,motif,haplotype,size,methylation,interruption_motif,read_support}
        if read_characteristics.motif.trim().is_empty() {
            continue;
        }
        
        let mut tmp_window_struct = per_window_struct {
            read_name: a.record().name(),
            window_start: window_start,
            window_end: window_end,
            motif: read_characteristics.motif,
            haplotype: 1,
            size: read_characteristics.size,
            min_methylation: read_characteristics.min_methylation,
            max_methylation: read_characteristics.max_methylation,
            avg_methylation: read_characteristics.avg_methylation,
            interruption_motif: "",
        }

        window_aggregate.push(tmp_window_struct);
    }
    //TODO::check per window map to see if there are any potential sequences that have more read support than the threshold
    //TODO::check per window map to see if there are any potential sequences that have larger count than the threshold
    //TODO:: if they do return in required format, else return with found_flag set to 0
    //FORMAT::if they do, send it back to main in the format potential_sequence maps to {found_flag,ref_start,ref_end,window_start,window_end,motif,haplotype,size,min methylation,max methylation,avg methylation,interruption_motif,read_support}
    return window_aggregate;
}