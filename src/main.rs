//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Sabiq Chaudhary (schaudhary@oicr.on.ca)
//---------------------------------------------------------
#[macro_use]
extern crate rust_htslib;
extern crate clap;

use rust_htslib::{bam, faidx, bam::Read, bam::record::CigarStringView};
use clap::{App, SubCommand, Arg, value_t};
use std::fs::File;
use std::path::Path;
use std::io::Write;

mod methylation_detection;

mod str_discovery;

mod str_sizing;

fn get_extension_from_filename(filename: &str) -> Option<&str> {    
    Path::new(filename)        
    .extension()        
    .and_then(OsStr::to_str)}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn main() {
    let matches = App::new("STRr10ToolKit")
        .version("0.1")
        .author("Sabiq Chaudhary <schaudhary@oicr.on.ca>")
        .about("Toolkit to discover, size and give methylation characteristics of STRs from BAMs generated with r10 and a modern basecaller which gives modified base tags")
        .arg(Arg::with_name("bam")
            .short("b")
            .long("bam")
            .takes_value(true)
            .help("the input bam file"))
        .arg(Arg::with_name("reference")
            .short("ref")
            .long("reference")
            .takes_value(true)
            .help("the reference genome"))
        .arg(Arg::with_name("output_file_name")
            .short("o")
            .long("output_file_name")
            .takes_value(true)
            .help("the output file name without extension"))
        .arg(Arg::with_name("output_directory")
            .long("output_directory")
            .takes_value(true)
            .help("the directory where all output files will be generated"))
        .arg(Arg::with_name("min_ins_size")
            .long("min_ins_size")
            .takes_value(true)
            .help("the minimum number of bases in the insert to call as an STR. Default Minimum is 50bp."))
        .arg(Arg::with_name("is_phased")
            .long("is_phased")
            .takes_value(true)
            .help("is phasing data present in the bam?"))
        .arg(Arg::with_name("min_read_support")
            .long("min_read_support")
            .takes_value(true)
            .help("the minimum number of reads that support the insert call. Minimum is set at 3 by default."))
        .arg(Arg::with_name("chromosomes")
            .long("chromosomes")
            .takes_value(true)
            .help("the chromosomes to be searched for STRs. Alternatively, you can provide a list of chromosomes in text file separated by a whitespace"))
        .arg(Arg::with_name("clean")
            .short("c")
            .long("clean")
            .takes_value(true)
            .help("cleanup any intermediate files"))
        .get_matches();

    let min_ins_size = value_t!(matches, "min_ins_size", u32).unwrap_or(50);
    let min_read_support = value_t!(matches, "min_read_support", u32).unwrap_or(3);
    let output_file_name = matches.value_of("output_file_name").unwrap_or("STRResults");
    let input_bam = matches.value_of("bam").unwrap();
    let reference_genome = matches.value_of("reference").unwrap()
    let output_directory = matches.value_of("output_directory").unwrap()
    let reads_phased = value_t!(matches.value_of("is_phased").unwrap())
    let chr_list_opt = matches.value_of("chromosomes").unwrap()
    let clean_flag = value_t!(matches.value_of("clean").unwrap())

    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);

    // set up header and faidx for pulling reference sequence
    let faidx = faidx::Reader::from_path(reference_genome).expect("Could not read reference genome:");
    let tid = header_view.tid(chromosome_name.as_bytes()).unwrap();
    let chromosome_length = header_view.target_len(tid).unwrap() as usize;
    let mut chromosome_sequence = faidx.fetch_seq_string(chromosome_name, 0, chromosome_length).unwrap();
    chromosome_sequence.make_ascii_uppercase();
    let chromosome_bytes = chromosome_sequence.as_bytes();

    if get_extension_from_filename(chr_list_opt) == "txt" {
        if let Ok(lines) = read_lines(chr_list_opt) {
            // Consumes the iterator, returns an (Optional) String
            for line in lines {
                chr_list = line
            }
        }
    }
    else {
        chr_list = chr_list_opt
    }

    for chr in chr_list {
        bam.fetch( chr );
        for p in bam.pileup() {
            let pileup = p.unwrap();
        }
    }

    //let regions_of_interest = methylation_detection::detect_methylation();

    let str_sites = str_discovery::detect_loci();

    //let str_size = str_sizing::detect_size();
}

