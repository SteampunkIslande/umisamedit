use std::env;
use bam::RecordWriter;
use bam::record::tags::StringType;
use bam::{ RecordReader, record::tags::TagValue };

use argparse_rs;
use argparse_rs::{ ArgParser, ArgType };

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
        "tag",
        None,
        't',
        true,
        "Tag to add at the end of each read name",
        ArgType::Option
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
    let tag_name: String = args.get("tag").unwrap();

    let mut reader = bam::BamReader::from_path(input_filename, 4).unwrap();
    let header = reader.header();

    let mut writer = bam::BamWriter::from_path(output_filename, header.clone()).unwrap();

    let mut record = bam::Record::new();
    loop {
        match reader.read_into(&mut record) {
            Ok(true) => {}
            Ok(false) => {
                break;
            }
            Err(e) => panic!("{}", e),
        }
        match
            record
                .tags()
                .get(
                    tag_name
                        .as_bytes()[..2]
                        .try_into()
                        .expect("Tag names must be exactly 2 bytes wide!")
                )
        {
            Some(TagValue::String(umi, StringType::String)) => {
                record.set_name([record.name(), b"_", umi].concat());
            }
            None => {}
            _ => {}
        }
        match writer.write(&record) {
            Ok(()) => {}
            Err(e) => {
                panic!("Error writing record:\n{}", e);
            }
        }
    }
}
