//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Sabiq Chaudhary (schaudhary@oicr.on.ca)
//---------------------------------------------------------

use rust_htslib::{bam, faidx, bam::Read, bam::record::CigarStringView}; 

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

/*function will take in the alignment, the read start, read end, the motif. It will go through the region specified on the 
read and read the corresponding methylation tags. It will average out the methylation over the region specified and give
the average methylation for the region. but it will also give the max and min methylation. if there are no tags,
it will return -1 or N/A.*/
pub fn detect_methylation -> u32 (start_pos: u32, end_pos: u32, motif: &String, record: &bam::Record) {
    println!("Inside methylation detection");
}