//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Sabiq Chaudhary (schaudhary@oicr.on.ca)
//---------------------------------------------------------

#[macro_use]
extern crate rust_htslib;
extern crate substring;

use rust_htslib::{bam, faidx, bam::Read, bam::record::CigarStringView};
use substring::Substring;
use std::cmp;
use std::collections::HashMap;

mod methylation_detection;
mod str_sizing;

struct per_read_struct {
    pub ref_start: u32,
    pub ref_end: u32,
    pub motif: &str,
    pub haplotype: u32,
    pub size: u32,
    pub methylation: u32,
    pub interruption_motif: &str,
}

struct per_window_struct {
    pub ref_start: u32,
    pub ref_end: u32,
    pub window_start: u32,
    pub window_end: u32,
    pub motif: &str,
    pub haplotype: u32,
    pub size: u32,
    pub min_methylation: u32,
    pub max_methylation: u32,
    pub avg_methylation: u32,
    pub interruption_motif: &str,
    pub read_support: u32,
}

pub struct all_charac_struct {
    pub found_flag: bool,
    pub ref_start: u32,
    pub ref_end: u32,
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

    for (key, value) in &*map {
        if value > max_val {
            max_key = key;
            max_val = value;
        }
    }

    return (max_key,max_val);
}

pub fn detect_loci(window_start: u32, window_end: u32, alignments: rust_htslib::bam::pileup::Alignments<'_>) {
    println!("Inside STR Discovery");
    //also pass a previous window ref_start and ref_end parameter. if present calculation aligns to the same region, skip this window and move to the next.

    for a in alignments {
        if a.record().seq().len() == 0 {
            continue;
        }

        //TODO::check bam record for any repeats here
        let cigar_string = a.cigar();
        if cigar_string.trim().is_empty() {
            continue;
        }

        let mut event_length = 0; //this is present to count the number of bases in each 'M','I','S', or other CIGAR events
        let mut read_counter = 0; //this is present to keep track of number of aligned bases till we hit the region of interest in the CIGAR string
        for letter in cigar_string {
            if letter.is_numeric() {
                event_length = 10*event_length + letter.to_digit(10).unwrap();
            }
            else {
                match letter {
                    'I' => {
                        if event_length > 50 {
                            let mut sequence_of_interest = a.record().seq().substring(read_counter,event_length); //slice it from (read_counter, event_length)
                            let mut potential_str_sequence = decompose_string(&sequence_of_interest);
                            if potential_str_sequence.trim().is_empty() {
                                continue;
                            }
                            read_counter += event_length;
                            event_length = 0;
                            //TODO::call sizing and methylation function if there is a potential str, to calc characteristics
                            //TODO::send sequence of interest as parameter to the sizing
                            //TODO::send record with read start and end positions to calculate methylation
                            //FORMAT::store in per read map which is potential_sequence maps to {ref_start,ref_end,motif,haplotype,size,methylation,interruption_motif}
                        }
                    },
                    'S' => {
                        if event_length > 50 {
                            let mut sequence_of_interest = a.record().seq().substring(read_counter,event_length); //slice it from (read_counter, event_length)
                            let mut potential_str_sequence = decompose_string(&sequence_of_interest);
                            if potential_str_sequence.trim().is_empty() {
                                continue;
                            }
                            read_counter += event_length;
                            event_length = 0;
                            //TODO::call sizing and methylation function if there is a potential str, to calc characteristics
                            //TODO::send sequence of interest as parameter to the sizing
                            //TODO::send record with read start and end positions to calculate methylation
                            //FORMAT::store in per read map which is potential_sequence maps to {ref_start,ref_end,motif,haplotype,size,methylation,interruption_motif}
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

        //TODO::aggregate all under per window map which is potential_sequence maps to {ref_start,ref_end,window_start,window_end,motif,haplotype,size,methylation,interruption_motif,read_support}
    }
    //TODO::check per window map to see if there are any potential sequences that have more read support than the threshold
    //TODO::check per window map to see if there are any potential sequences that have larger count than the threshold
    //TODO:: if they do return in required format, else return with found_flag set to 0
    //FORMAT::if they do, send it back to main in the format potential_sequence maps to {found_flag,ref_start,ref_end,window_start,window_end,motif,haplotype,size,min methylation,max methylation,avg methylation,interruption_motif,read_support}
}