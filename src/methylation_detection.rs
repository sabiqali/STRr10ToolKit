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
pub fn detect_methylation -> u32 (read_start: u32, read_end: u32, record: &bam::Record) {
    //TODO::use the MM and ML tags in the region from the above coordinates to calculate methylation for this record. return it. 
    //ADDITIONAL::Also, get the haplotype to get which haplotype the read belongs to. 

    if let Ok(Aux::String(mm_str)) = get_mm_tag(record) {
        if let Ok(Aux::ArrayU8(probability_array)) = get_ml_tag(record) {

        }
    }
}