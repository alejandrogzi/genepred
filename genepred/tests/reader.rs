use genepred::reader::Reader;
use genepred::{Bed12, Bed3, Bed4, Bed6};

#[test]
fn test_reader_from_string_bed3() {
    let data = "chr1\t10\t20\nchr1\t30\t40";
    let mut reader: Reader<Bed3> =
        Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 2);

    let first = &records[0];
    assert_eq!(first.chrom(), b"chr1".as_ref());
    assert_eq!(first.start(), 10);
    assert_eq!(first.end(), 20);

    let second = &records[1];
    assert_eq!(second.chrom(), b"chr1".as_ref());
    assert_eq!(second.start(), 30);
    assert_eq!(second.end(), 40);
}

#[test]
fn test_reader_from_string_bed4() {
    let data = "chr1\t10\t20\tgeneA\nchr1\t30\t40\tgeneB";
    let mut reader: Reader<Bed4> =
        Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 2);

    let first = &records[0];
    assert_eq!(first.chrom(), b"chr1".as_ref());
    assert_eq!(first.start(), 10);
    assert_eq!(first.end(), 20);
    assert_eq!(first.name().unwrap(), b"geneA".as_ref());

    let second = &records[1];
    assert_eq!(second.chrom(), b"chr1".as_ref());
    assert_eq!(second.start(), 30);
    assert_eq!(second.end(), 40);
    assert_eq!(second.name().unwrap(), b"geneB".as_ref());
}

#[test]
fn test_reader_from_string_bed6() {
    let data = "chr1\t10\t20\tgeneA\t100\t+\nchr1\t30\t40\tgeneB\t200\t-";
    let mut reader: Reader<Bed6> =
        Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 2);

    let first = &records[0];
    assert_eq!(first.chrom(), b"chr1".as_ref());
    assert_eq!(first.start(), 10);
    assert_eq!(first.end(), 20);
    assert_eq!(first.name().unwrap(), b"geneA".as_ref());
    assert_eq!(first.score().unwrap(), 100 as u16);
    assert_eq!(first.strand().unwrap().to_string(), "+");

    let second = &records[1];
    assert_eq!(second.chrom(), b"chr1".as_ref());
    assert_eq!(second.start(), 30);
    assert_eq!(second.end(), 40);
    assert_eq!(second.name().unwrap(), b"geneB".as_ref());
    assert_eq!(second.score().unwrap(), 200 as u16);
    assert_eq!(second.strand().unwrap().to_string(), "-");
}

#[test]
fn test_reader_from_string_bed12() {
    let data = "chr1\t10\t100\tgeneA\t1000\t+\t10\t100\t255,0,0\t2\t10,20,\t0,30,\n";
    let mut reader: Reader<Bed12> =
        Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 1);

    let first = &records[0];
    assert_eq!(first.chrom(), b"chr1".as_ref());
    assert_eq!(first.start(), 10);
    assert_eq!(first.end(), 100);
    assert_eq!(first.name().unwrap(), b"geneA".as_ref());
    assert_eq!(first.score().unwrap(), 1000 as u16);
    assert_eq!(first.strand().unwrap().to_string(), "+");
    assert_eq!(first.thick_start().unwrap(), 10);
    assert_eq!(first.thick_end().unwrap(), 100);
    assert_eq!(first.item_rgb().unwrap().to_string(), "255,0,0");
    assert_eq!(first.block_count().unwrap(), 2);
    assert_eq!(first.block_sizes().unwrap(), vec![10, 20]);
    assert_eq!(first.block_starts().unwrap(), vec![0, 30]);
}

#[test]
fn test_reader_invalid_line() {
    let data = "chr1\t10\t20\nmalformed_line\nchr2\t50\t60";
    let mut reader: Reader<Bed3> =
        Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
    let records: Vec<_> = reader.records().collect();
    assert_eq!(records.len(), 3);

    assert!(records[0].is_ok());
    assert!(records[1].is_err());
    assert!(records[2].is_ok());
}

#[test]
fn test_reader_empty_input() {
    let data = "";
    let mut reader: Reader<Bed3> =
        Reader::from_reader(std::io::Cursor::new(data.as_bytes())).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert!(records.is_empty());
}
