use genepred::reader::Reader;
use genepred::{Bed12, Bed3, Bed4, Bed6, ExtraValue, Gff, Gtf, Strand};

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
    assert_eq!(first.strand().unwrap().to_string(), "+");

    let second = &records[1];
    assert_eq!(second.chrom(), b"chr1".as_ref());
    assert_eq!(second.start(), 30);
    assert_eq!(second.end(), 40);
    assert_eq!(second.name().unwrap(), b"geneB".as_ref());
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
    assert_eq!(first.strand().unwrap().to_string(), "+");
    assert_eq!(first.thick_start().unwrap(), 10);
    assert_eq!(first.thick_end().unwrap(), 100);
    assert_eq!(first.block_count().unwrap(), 2);
    assert_eq!(first.block_starts().unwrap(), vec![10, 40]);
    assert_eq!(first.block_ends().unwrap(), vec![20, 60]);
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

#[test]
fn test_reader_gxf_from_path() {
    let path = "tests/data/simple.gtf";
    let mut reader: Reader<Gtf> = Reader::from_path(path).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();

    assert_eq!(records.len(), 1);
    let gene = &records[0];
    assert_eq!(gene.chrom(), b"chr1".as_ref());
    assert_eq!(gene.start(), 99);
    assert_eq!(gene.end(), 200);
    assert_eq!(gene.name().unwrap(), b"GeneOne".as_ref());
    assert_eq!(gene.strand().unwrap(), Strand::Forward);
    assert_eq!(gene.block_count().unwrap(), 2);
    assert_eq!(gene.block_starts().unwrap(), &[99, 169]);
    assert_eq!(gene.block_ends().unwrap(), &[150, 200]);
    assert_eq!(gene.thick_start().unwrap(), 119);
    assert_eq!(gene.thick_end().unwrap(), 180);
}

#[test]
fn test_reader_bed3_from_path() {
    let path = "tests/data/bed3.bed";
    let mut reader: Reader<Bed3> = Reader::from_path(path).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].start(), 0);
    assert_eq!(records[0].end(), 100);
    assert_eq!(records[1].start(), 150);
    assert_eq!(records[1].end(), 200);
}

#[test]
fn test_reader_bed6_from_path() {
    let path = "tests/data/bed6.bed";
    let mut reader: Reader<Bed6> = Reader::from_path(path).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 2);

    let first = &records[0];
    assert_eq!(first.name().unwrap(), b"geneA".as_ref());
    assert_eq!(first.strand().unwrap(), Strand::Forward);

    let second = &records[1];
    assert_eq!(second.name().unwrap(), b"geneB".as_ref());
    assert_eq!(second.strand().unwrap(), Strand::Reverse);
}

#[test]
fn test_reader_bed12_from_path() {
    let path = "tests/data/bed12.bed";
    let mut reader: Reader<Bed12> = Reader::from_path(path).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();

    assert_eq!(records.len(), 1);
    let gene = &records[0];
    assert_eq!(gene.name().unwrap(), b"txA".as_ref());
    assert_eq!(gene.block_count().unwrap(), 2);
    assert_eq!(gene.block_starts().unwrap(), &[100, 300]);
    assert_eq!(gene.block_ends().unwrap(), &[180, 360]);
    assert_eq!(gene.thick_start().unwrap(), 120);
    assert_eq!(gene.thick_end().unwrap(), 360);
}

#[test]
fn test_reader_bed12_with_additional_fields() {
    let path = "tests/data/bed12_extra.bed";
    let mut reader: Reader<Bed12> = Reader::from_path_with_additional_fields(path, 2).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();

    assert_eq!(records.len(), 1);
    let gene = &records[0];
    assert_eq!(gene.name().unwrap(), b"txB".as_ref());
    assert_eq!(gene.block_count().unwrap(), 3);
    assert_eq!(gene.block_starts().unwrap(), &[50, 150, 200]);
    assert_eq!(gene.block_ends().unwrap(), &[110, 190, 250]);

    let extras = gene.extras();
    match extras.get(&b"13".to_vec()) {
        Some(ExtraValue::Scalar(value)) => assert_eq!(value, b"foo"),
        other => panic!("unexpected extras[13]: {:?}", other),
    }
    match extras.get(&b"14".to_vec()) {
        Some(ExtraValue::Scalar(value)) => assert_eq!(value, b"bar"),
        other => panic!("unexpected extras[14]: {:?}", other),
    }
}

#[test]
fn test_reader_gff_from_path() {
    let path = "tests/data/simple.gff";
    let mut reader: Reader<Gff> = Reader::from_path(path).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();

    assert_eq!(records.len(), 1);
    let gene = &records[0];
    assert_eq!(gene.name().unwrap(), b"Tx1".as_ref());
    assert_eq!(gene.chrom(), b"chr1".as_ref());
    assert_eq!(gene.start(), 99);
    assert_eq!(gene.end(), 200);
    assert_eq!(gene.strand().unwrap(), Strand::Forward);
    assert_eq!(gene.block_count().unwrap(), 2);
    assert_eq!(gene.block_starts().unwrap(), &[99, 169]);
    assert_eq!(gene.block_ends().unwrap(), &[150, 200]);
    assert_eq!(gene.thick_start().unwrap(), 119);
    assert_eq!(gene.thick_end().unwrap(), 180);
}

#[cfg(feature = "compression")]
#[test]
fn test_reader_bed3_gz_from_path() {
    let path = "tests/data/bed3.bed.gz";
    let mut reader: Reader<Bed3> = Reader::from_path(path).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].start(), 0);
    assert_eq!(records[1].end(), 200);
}

#[cfg(feature = "compression")]
#[test]
fn test_reader_gtf_gz_from_path() {
    let path = "tests/data/simple.gtf.gz";
    let mut reader: Reader<Gtf> = Reader::from_path(path).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 1);
    let gene = &records[0];
    assert_eq!(gene.name().unwrap(), b"GeneOne".as_ref());
    assert_eq!(gene.block_count().unwrap(), 2);
}
