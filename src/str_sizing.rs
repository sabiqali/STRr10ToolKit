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

/*function will take in the alignment, the read start, read end, the motif. it will search through the insert positions
on the read for the motif in a simplistic approach, as we assume the accuracy from r10 kit 14 onwards is high. 
it will return a dict with a count of the supplied motif and if there are any interruptions in the motif.*/
pub fn detect_size -> u32 (start_pos: u32, end_pos: u32, motif: &String, record: &bam::Record) {

    println!("Inside STR sizing");
}