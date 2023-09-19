use std::collections::HashMap;
use std::rc::Rc;
use linear_map::LinearMap;
use rust_htslib::bam::HeaderView;
use rust_htslib::bam::header::HeaderRecord;
use std::env;
use std::error::Error;

use rust_htslib::bam::record::Aux;
use rust_htslib::{ bam, bam::Read, bam::Record };

use argparse_rs;
use argparse_rs::{ ArgParser, ArgType };
use csv;

use serde::{ Deserialize, Serialize };

#[derive(Serialize, Deserialize, Debug)]
struct CsvRecord {
    key: String,
    value: String,
}

fn read_csv_to_hashmap(filename: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
    // Create a CSV reader
    let mut rdr = csv::ReaderBuilder
        ::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(filename)?;

    // Create a HashMap to store the key-value pairs
    let mut hashmap = HashMap::new();

    // Iterate through each record in the CSV
    for result in rdr.deserialize() {
        let record: CsvRecord = result?;
        // Insert the key-value pair into the HashMap
        hashmap.insert(record.key, record.value);
    }

    Ok(hashmap)
}

fn rename_contig(
    header: &mut bam::Header,
    original_header_map: &HashMap<String, Vec<LinearMap<String, String>>>,
    translation_map: HashMap<String, String>
) {
    for (key, val) in original_header_map {
        for names in val {
            let mut new_record = HeaderRecord::new(key.as_bytes());
            //If tag_value doesn't have a translation, don't add the tag
            let mut ignore_tag = false;
            for name in names {
                let (tag_name, mut tag_value) = name;

                if key == "SQ" && tag_name == "SN" {
                    tag_value = match translation_map.get(tag_value) {
                        Some(translation) => { translation }
                        None => {
                            ignore_tag = true;
                            tag_value
                        }
                    };
                }
                if !ignore_tag {
                    new_record.push_tag(tag_name.as_bytes(), tag_value);
                }
            }
            if !ignore_tag {
                header.push_record(&new_record);
            }
        }
    }

    // Below, code to make sure the header was edited correctly
    // for (key, val) in &header.to_hashmap() {
    //     for names in val {
    //         print!("@{}\t", key);
    //         for name in names {
    //             print!("{}:{}\t", name.0, name.1);
    //         }
    //         print!("\n");
    //     }
    // }
}

fn main() {
    let raw_args: Vec<String> = env::args().collect();

    let mut parser = ArgParser::new("argparse".into());

    parser.add_opt("input", None, 'i', true, "Input file to edit.", ArgType::Positional(0));
    parser.add_opt(
        "output",
        None,
        'o',
        true,
        "Output file to write edited bam to.",
        ArgType::Positional(1)
    );
    parser.add_opt(
        "translation",
        None,
        't',
        true,
        "Translation file to rename contigs.",
        ArgType::Positional(2)
    );

    let args = match parser.parse(raw_args.iter()) {
        Ok(args) => { args }
        Err(e) => {
            parser.help();
            panic!("Error parsing arguments:\n{e}\n");
        }
    };

    let input_filename: String = args.get("input").unwrap();
    let output_filename: String = args.get("output").unwrap();
    let translation_filename: String = args.get("translation").unwrap();

    let tag_name: String = String::from("RX");

    let mut reader = bam::Reader::from_path(&input_filename).unwrap();

    let original_header = bam::Header::from_template(reader.header());

    let mut new_header = bam::Header::new();

    let translation_map = read_csv_to_hashmap(&translation_filename).expect(
        "Translation file should exist and have exactly two columns"
    );

    //Creates a new header, renaming the contigs using translation_map
    rename_contig(&mut new_header, &original_header.to_hashmap(), translation_map);

    let mut writer = bam::Writer
        ::from_path(output_filename, &new_header, bam::Format::Bam)
        .unwrap();

    let mut record = Record::new();
    let rc_new_headerview = Rc::new(HeaderView::from_header(&new_header));

    while let Some(_) = reader.read(&mut record) {
        match record.aux(tag_name.as_bytes()) {
            Ok(value) => {
                if let Aux::String(v) = value {
                    record.set_qname(&[record.qname(), b"_", v.as_bytes()].concat());
                }
            }
            Err(..) => {}
        }

        if record.tid() < rc_new_headerview.target_count().try_into().unwrap() {
            record.set_header(rc_new_headerview.clone());
            if record.tid() >= rc_new_headerview.target_count().try_into().unwrap() {
                println!("Sérieux?");
            }
            match writer.write(&record) {
                Ok(()) => {}
                Err(e) => {
                    println!("{}", e);
                }
            }
        } else {
            println!("We have an imposter {:?}", record);
            break;
        }
    }
}
