//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Sabiq Chaudhary (schaudhary@oicr.on.ca)
//---------------------------------------------------------

#[macro_use]
extern crate rust_htslib;

use rust_htslib::{bam, faidx, bam::Read, bam::record::CigarStringView, bam::record::Aux}; 

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
pub fn detect_methylation -> (f64,f64,f64) (read_start: u32, read_end: u32, record: &bam::Record) {
    //TODO::use the MM and ML tags in the region from the above coordinates to calculate methylation for this record. return it. 
    //ADDITIONAL::Also, get the haplotype to get which haplotype the read belongs to. 

    let mut read_tag_count = 0;
    let mut read_tag_current_base = 0;

    let mut max_methylation = 0;
    let mut min_methylation = 0;
    let mut avg_methylation = 0;
    let mut total_methylation = 0;

    //does this need to be nested?
    Aux::String(mm_str)) = get_mm_tag(record) 
    Aux::ArrayU8(probability_array)) = get_ml_tag(record) 
            
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