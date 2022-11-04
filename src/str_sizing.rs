//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Sabiq Chaudhary (schaudhary@oicr.on.ca)
//---------------------------------------------------------

#[macro_use]
extern crate rust_htslib;
extern crate substring;

use rust_htslib::{bam, faidx, bam::Read, bam::record::CigarStringView};
use substring::Substring;

pub struct size_struct {
    pub count: u32,
    pub interruption_motif: String,
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

/*function will take in the alignment, the read start, read end, the motif. it will search through the insert positions
on the read for the motif in a simplistic approach, as we assume the accuracy from r10 kit 14 onwards is high. 
it will return a dict with a count of the supplied motif and if there are any interruptions in the motif.*/
pub fn detect_size -> size_struct (sequence_of_interest: &str, potential_str_sequence: &str) {
    //TODO::this function will go throught the sequence of interest and check the count and return it
    //additional:: check if there are any interruptions while counting. if we do this, report count and interruption.

    let v: Vec<_> = sequence_of_interest.match_indices(potential_str_sequence).collect();
    if v.len() <= 2 {
        return 0;
    }

    let mut motif_length = potential_str_sequence.len()

    let mut output_variable = size_struct;

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
        i++;
        j++;
    }
    return output_variable;
}