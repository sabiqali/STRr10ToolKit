//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Sabiq Chaudhary (schaudhary@oicr.on.ca)
//---------------------------------------------------------

extern crate rust_htslib;
extern crate substring;
//use crate::str_sizing::size_struct;

use rust_htslib::{bam, faidx, bam::Read, bam::record::CigarStringView, bam::record::Aux};
use substring::Substring;
use std::cmp;
use std::collections::HashMap;
use std::io::BufReader;

//mod methylation_detection;
//mod str_sizing;

struct per_read_struct {
    pub read_start: u32,
    pub read_end: u32,
    pub motif: String,
    pub haplotype: u32,
    pub size: u32,
    pub min_methylation: f64,
    pub max_methylation: f64,
    pub avg_methylation: f64,
    pub interruption_motif: String,
}

pub struct per_window_struct {
    pub read_name: String,
    pub window_start: u32,
    pub window_end: u32,
    pub motif: String,
    pub haplotype: u32,
    pub size: u32,
    pub min_methylation: f64,
    pub max_methylation: f64,
    pub avg_methylation: f64,
    pub interruption_motif: String,
}

pub struct all_charac_struct {
    pub found_flag: bool,
    pub window_start: u32,
    pub window_end: u32,
    pub motif: String,
    pub haplotype: u32,
    pub size: u32,
    pub min_methylation: f64,
    pub max_methylation: f64,
    pub avg_methylation: f64,
    pub interruption_motif: String,
    pub read_support: u32,
}

pub struct decomposer_struct {
    pub potential_sequences_in_window: Vec<String>,
    pub potential_count_in_window: Vec<u32>,
}

/// Returns the length of the longest substring without repeating characters.
///
/// # Arguments
///
/// * 's' - The input string
/*pub fn length_of_longest_substring(s: String) -> i32 {
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
}*/

pub struct size_struct {
    pub count: u32,
    pub interruption_motif: String,
}

/*function will take in the alignment, the read start, read end, the motif. it will search through the insert positions
on the read for the motif in a simplistic approach, as we assume the accuracy from r10 kit 14 onwards is high. 
it will return a dict with a count of the supplied motif and if there are any interruptions in the motif.*/
pub fn detect_size(sequence_of_interest: &str, potential_str_sequence: &str) -> size_struct  {
    //TODO::this function will go throught the sequence of interest and check the count and return it
    //additional:: check if there are any interruptions while counting. if we do this, report count and interruption.

    let v: Vec<_> = sequence_of_interest.match_indices(potential_str_sequence).collect();
    if v.len() <= 2 {
        return size_struct {
            count: 0,
            interruption_motif: "".to_string(),
        };
    }

    let mut motif_length = potential_str_sequence.len();

    let mut output_variable: size_struct;

    let mut i=0;
    let mut j=1;
    let mut count = 0;
    while j <= v.len() {
        if v[j] - v[i] == motif_length {
            output_variable.count+=1;
        }
        else {
            output_variable.interruption_motif = sequence_of_interest.substring(v[i], v[i]+motif_length)
        }
        i+=1;
        j+=1;
    }
    return output_variable;
}

fn get_haplotag_from_record(record: &bam::Record) -> Option<i32> {
    match record.aux(b"HP") {
        Ok(value) => {
            if let Aux::I32(v) = value {
                return Some(v)
            } else if let Aux::U8(v) = value {
                return Some(v as i32) // whatshap encodes with U8
            } else {
                return None
            }
        }
        Err(_e) => return None
    }
}

pub fn get_mm_tag(record: &bam::Record) -> Result<Aux, rust_htslib::errors::Error> {
    let r = match record.aux(b"MM") {
        Ok(v) => Ok(v),
        Err(_) => record.aux(b"Mm")
    };
    return r;
}

pub fn get_ml_tag(record: &bam::Record) -> Result<Aux, rust_htslib::errors::Error> {
    let r = match record.aux(b"ML") {
        Ok(v) => Ok(v),
        Err(_) => record.aux(b"Ml")
    };
    return r;
}

/*function will take in the alignment, the read start, read end, the motif. It will go through the region specified on the 
read and read the corresponding methylation tags. It will average out the methylation over the region specified and give
the average methylation for the region. but it will also give the max and min methylation. if there are no tags,
it will return -1 or N/A.*/
pub fn detect_methylation(read_start: u32, read_end: u32, record: &bam::Record) -> (f64,f64,f64){
    //TODO::use the MM and ML tags in the region from the above coordinates to calculate methylation for this record. return it. 
    //ADDITIONAL::Also, get the haplotype to get which haplotype the read belongs to. 

    let mut read_tag_count = 0;
    let mut read_tag_current_base = 0;

    let mut max_methylation = 0;
    let mut min_methylation = 0;
    let mut avg_methylation = 0;
    let mut total_methylation = 0;

    //does this need to be nested?
    let mut tmp_mm_str = get_mm_tag(record); 
    let mut tmp_probability_array = get_ml_tag(record); 

    let mut mm_str = Aux::String(tmp_mm_str);
    let mut probability_array = Aux::ArrayU8(tmp_probability_array); 
            
    // TODO: handle multiple mods
    if mm_str.matches(';').count() != 1 {
        return None;
    }

    let first_mod_str = mm_str.split(';').next().unwrap();
    let mod_meta = first_mod_str.split(',').next().unwrap().as_bytes();

    //derived from mbtools(https://github.com/jts/mbtools/blob/b61e0c15b2d058d3b6b55fc3f46d6d64eb674d15/src/main.rs#L213)
    //Use own version of parsing by comma if it doesn't work.
    for (token, encoded_probability) in first_mod_str.split(',').skip(1).zip(probability_array.iter()) {
        let mut position = token.parse::<usize>().unwrap();
        let methylation_base_probability = encoded_probability as f64 / 255.0;

        read_tag_count = read_tag_count + position + 1;
        if read_tag_count >= read_start && read_tag_count <= read_end {
            if methylation_base_probability > max_methylation {
                max_methylation = methylation_base_probability;
            }
            if methylation_base_probability < min_methylation {
                min_methylation = methylation_base_probability;
            }
            total_methylation += methylation_base_probability;
        }
    }
    if read_tag_count > read_end {
        avg_methylation = total_methylation / (read_end - read_start);
    }
    else if read_tag_count < read_end && read_tag_count > read_start {
        avg_methylation = total_methylation / (read_tag_count - read_start);
    }
    else if read_tag_count < read_start {
        avg_methylation = 0
    }
    return (avg_methylation,min_methylation, max_methylation);
}

pub fn decompose_string(sequence_of_interest: &str, lower_limit: u32, upper_limit: u32) -> (&str,u32) {
    //TODO::take in the string of interest from the CIGAR String and then find if there is a repeating unit present

    let mut decomposer_return_struct: decomposer_struct;

    let mut subsequences: HashMap<&str, u32> = HashMap::new();

    for motif_length in (lower_limit,upper_limit) {

        let mut lower_window_var = 0;
        let mut upper_window_var = motif_length;

        while upper_window_var <= sequence_of_interest.len() {
            let mut sequence_in_window = sequence_of_interest.substring(lower_window_var,upper_window_var);

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

pub fn detect_loci(window_start: u32, window_end: u32, alignments: rust_htslib::bam::pileup::Alignments<'_>, lower_limit: u32, upper_limit: u32, min_read_support: u32, min_ins_size: u32, discovery_sensitivity: u32) -> Vec<per_window_struct> {
    println!("Inside STR Discovery");
    //also pass a previous window ref_start and ref_end parameter. if present calculation aligns to the same region, skip this window and move to the next.

    let mut window_aggregate: Vec<per_window_struct> = Vec::new();

    //let mut read_support_counter: u32;

    for a in alignments {
        if a.record().seq().len() == 0 {
            continue;
        }

        //TODO::check bam record for any repeats here
        let mut cigar_string = a.record().cigar();
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
                            let (mut potential_str_sequence,mut potential_count_from_discovery) = decompose_string(&sequence_of_interest,lower_limit,upper_limit);
                            if potential_str_sequence.trim().is_empty() {
                                continue;
                            }
                            //TODO::call sizing and methylation function if there is a potential str, to calc characteristics
                            //TODO::send sequence of interest as parameter to the sizing
                            //TODO::send record with read start and end positions to calculate methylation
                            //FORMAT::store in per read map which is potential_sequence maps to {ref_start,ref_end,motif,haplotype,size,methylation,interruption_motif}
                            
                            let mut size_func_return: size_struct;
                            size_func_return = detect_size(&sequence_of_interest,&potential_str_sequence);

                            if size_func_return.count <= discovery_sensitivity {
                                continue;
                            }

                            let (avg_methylation,min_methylation, max_methylation) = detect_methylation(read_counter, read_counter+event_length, &a.record());

                            read_characteristics.read_start = read_counter;
                            read_characteristics.read_end = read_counter+event_length;
                            read_characteristics.motif = potential_str_sequence.to_string();
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
                            let (mut potential_str_sequence,mut potential_count_from_discovery) = decompose_string(&sequence_of_interest,lower_limit,upper_limit);
                            if potential_str_sequence.trim().is_empty() {
                                continue;
                            }
                            //TODO::call sizing and methylation function if there is a potential str, to calc characteristics
                            //TODO::send sequence of interest as parameter to the sizing
                            //TODO::send record with read start and end positions to calculate methylation
                            //FORMAT::store in per read map which is read_name maps to {potential_motif,ref_start,ref_end,motif,haplotype,size,methylation,interruption_motif}
                            let mut size_func_return: size_struct;
                            size_func_return = detect_size(&sequence_of_interest,&potential_str_sequence);

                            if size_func_return.count <= discovery_sensitivity {
                                continue;
                            }

                            let (avg_methylation,min_methylation, max_methylation) = detect_methylation(read_counter, read_counter+event_length, &a.record());
                            
                            read_characteristics.read_start = read_counter;
                            read_characteristics.read_end = read_counter+event_length;
                            read_characteristics.motif = potential_str_sequence.to_string();
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
            read_name: a.record().qname(),
            window_start: window_start,
            window_end: window_end,
            motif: read_characteristics.motif,
            haplotype: 1,
            size: read_characteristics.size,
            min_methylation: read_characteristics.min_methylation,
            max_methylation: read_characteristics.max_methylation,
            avg_methylation: read_characteristics.avg_methylation,
            interruption_motif: "".to_string(),
        };

        window_aggregate.push(tmp_window_struct);
    }
    //TODO::check per window map to see if there are any potential sequences that have more read support than the threshold
    //TODO::check per window map to see if there are any potential sequences that have larger count than the threshold
    //TODO:: if they do return in required format, else return with found_flag set to 0
    //FORMAT::if they do, send it back to main in the format potential_sequence maps to {found_flag,ref_start,ref_end,window_start,window_end,motif,haplotype,size,min methylation,max methylation,avg methylation,interruption_motif,read_support}
    return window_aggregate;
}