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

use std::str;

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

fn header_rename_contig(
    header: &mut bam::Header,
    new_tid_map: &mut HashMap<i32, i32>,
    original_header_map: &HashMap<String, Vec<LinearMap<String, String>>>,
    translation_map: HashMap<String, String>
) {
    assert!(new_tid_map.is_empty(), "New TID map should start empty!");
    for (key, val) in original_header_map {
        //Keep track of the mapping of kept contigs vs skipped contigs, important when writing records
        let mut contig_counter = 0;
        let mut kept_contigs_counter = 0;

        for names in val {
            let mut new_record = HeaderRecord::new(key.as_bytes());
            //If tag_value doesn't have a translation, don't add the tag
            let mut ignore_tag = false;

            for name in names {
                let (tag_name, mut tag_value) = name;

                if key == "SQ" && tag_name == "SN" {
                    tag_value = match translation_map.get(tag_value) {
                        Some(translation) => {
                            match new_tid_map.insert(contig_counter, kept_contigs_counter) {
                                Some(_) => {
                                    panic!(
                                        "Logical error: how did we end up returning back to a previous value of an incremental?"
                                    );
                                }
                                None => {}
                            }
                            kept_contigs_counter += 1;
                            translation
                        }
                        None => {
                            ignore_tag = true;
                            tag_value
                        }
                    };
                    contig_counter += 1;
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

    // UNCOMMENT TO DEBUG
    // println!("{:?}", new_tid_map);
    // let hv: HeaderView = HeaderView::from_header(header);
    // for (k, v) in new_tid_map {
    //     println!(
    //         "Before: {}->{}\tAfter:{}->{:#?}",
    //         k,
    //         original_header_map["SQ"]
    //             .iter()
    //             .nth(*k as usize)
    //             .unwrap()["SN"],
    //         v,
    //         String::from_utf8(hv.tid2name(*v as u32).to_vec()).unwrap()
    //     );
    // }
    // panic!();

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
    let mut new_tid_map = HashMap::<i32, i32>::new();

    //Creates a new header, renaming the contigs using translation_map
    header_rename_contig(
        &mut new_header,
        &mut new_tid_map,
        &original_header.to_hashmap(),
        translation_map
    );

    let mut writer = bam::Writer
        ::from_path(output_filename, &new_header, bam::Format::Bam)
        .unwrap();

    let mut record = Record::new();
    let rc_new_headerview = Rc::new(HeaderView::from_header(&new_header));

    // Just two variables to get information about skipped contigs
    let mut skipped_contig_names: HashMap<String, u32> = HashMap::new();
    let mut skipped_records: u32 = 0;

    while let Some(_) = reader.read(&mut record) {
        match record.aux(tag_name.as_bytes()) {
            Ok(value) => {
                if let Aux::String(v) = value {
                    record.set_qname(&[record.qname(), b"_", v.as_bytes()].concat());
                }
            }
            Err(..) => {}
        }

        //Try to get the mapped contig ID
        match (new_tid_map.get(&record.tid()), new_tid_map.get(&record.mtid())) {
            //This contig ID was kept during translation
            (Some(new_tid), Some(new_mtid)) => {
                record.set_header(rc_new_headerview.clone());
                record.set_tid(*new_tid);
                record.set_mtid(*new_mtid);
                match writer.write(&record) {
                    Ok(()) => {}
                    Err(e) => {
                        eprintln!("{}", e);
                    }
                }
            }
            (Some(new_tid), None) => {
                record.set_header(rc_new_headerview.clone());
                record.set_tid(*new_tid);
                record.set_mtid(-1);
                match writer.write(&record) {
                    Ok(()) => {}
                    Err(e) => {
                        eprintln!("{}", e);
                    }
                }
            }
            (None, Some(new_mtid)) => {
                record.set_header(rc_new_headerview.clone());
                record.set_tid(-1);
                record.set_mtid(*new_mtid);
                match writer.write(&record) {
                    Ok(()) => {}
                    Err(e) => {
                        eprintln!("{}", e);
                    }
                }
            }
            //This contig ID was not kept in translation
            (None, None) => {
                let skipped_contig_name = String::from_utf8(
                    reader.header().tid2name(record.tid().try_into().unwrap()).to_vec()
                ).unwrap();
                let skipped_contig_count = skipped_contig_names.get(&skipped_contig_name);
                match skipped_contig_count {
                    Some(val) => {
                        skipped_contig_names.insert(skipped_contig_name, val + 1);
                    }
                    None => {
                        skipped_contig_names.insert(skipped_contig_name, 1);
                    }
                }
                skipped_records += 1;
            }
        }
    }
    eprintln!("Done properly renaming contigs, skipped {} records.", skipped_records);
    eprintln!("Skipped contigs {:?}", skipped_contig_names);
}
