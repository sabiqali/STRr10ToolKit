//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Sabiq Chaudhary (schaudhary@oicr.on.ca)
//---------------------------------------------------------
//#[macro_use]
extern crate rust_htslib;

use rust_htslib::{bam, faidx, bam::Read, bam::record::CigarStringView};
use clap::{App, SubCommand, Arg, value_t};
use itertools::Itertools;
use std::fs::File;
use std::path::Path;
use std::io::Write;

mod methylation_detection;

mod str_discovery;

mod str_sizing;

fn main() {
    let matches = App::new("STRr10ToolKit")
        .version("0.1")
        .author("Sabiq Chaudhary <schaudhary@oicr.on.ca>")
        .about("Toolkit to discover, size and give methylation characteristics from BAMs generated with r10 and a modern basecaller which gives modified base tags")
        .arg(Arg::with_name("bam")
            .short("b")
            .long("bam")
            .takes_value(true)
            .help("the input bam file"))
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
        .arg(Arg::with_name("clean")
            .short("c")
            .long("clean")
            .takes_value(true)
            .help("cleanup any intermediate files"));

    let mut min_ins_size = value_t!(matches, "min_ins_size", u32).unwrap_or(50);

    let mut output_file = matches.value_of("output_file_name").unwrap_or("STRResults");

    println!("{} {}",min_ins_size,output_file);
    //let mut 

    //let mut regions_of_interest = methylation_detection::detect_STRs()
}
