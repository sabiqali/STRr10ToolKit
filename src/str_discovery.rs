//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Sabiq Chaudhary (schaudhary@oicr.on.ca)
//---------------------------------------------------------

mod methylation_detection;
mod str_sizing;

pub fn decompose_string(sequence_of_interest: &str) {
    //TODO::take in the string of interest from the CIGAR String and then find if there is a repeating unit present
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
                            let mut sequence_of_interest = a.record().seq(); //slice it from (read_counter, event_length)
                            let mut potential_str_sequence = decompose_string(&sequence_of_interest);
                            if potential_str_sequence.trim().is_empty() {
                                continue;
                            }
                            //TODO::call sizing and methylation function if there is a potential str, to calc characteristics
                            //TODO::send sequence of interest as parameter to the sizing
                            //TODO::send record with read start and end positions to calculate methylation
                            //FORMAT::store in per read map which is potential_sequence maps to {ref_start,ref_end,motif,haplotype,size,methylation}
                        }
                    },
                    'S' => {

                    },
                    _ => {
                        read_counter += event_length;
                        event_length = 0;
                        continue;    
                    },
                }
            }
        }

        //TODO::aggregate all under per window map which is potential_sequence maps to {ref_start,ref_end,window_start,window_end,motif,haplotype,size,methylation,read_support}
    }
    //TODO::check per window map to see if there are any potential sequences that have more read support than the threshold
    //TODO::check per window map to see if there are any potential sequences that have larger count than the threshold
    //TODO:: if they do return in required format, else return with found_flag set to 0
    //FORMAT::if they do, send it back to main in the format potential_sequence maps to {found_flag,ref_start,ref_end,window_start,window_end,motif,haplotype,size,min methylation,max methylation,avg methylation,read_support}
}