use genepred::{
    genepred::{ExtraValue, Extras, GenePred},
    strand::Strand,
    Bed12, Bed3, Gff, Gtf, Reader, Writer, WriterOptions,
};

#[test]
fn write_gtf_from_genepred() {
    let mut extras = Extras::new();
    extras.insert(
        b"tag".to_vec(),
        ExtraValue::Array(vec![b"a".to_vec(), b"b".to_vec()]),
    );
    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 99, 200, extras);
    gene.set_name(Some(b"tx1".to_vec()));
    gene.set_strand(Some(Strand::Forward));
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![99, 169]));
    gene.set_block_ends(Some(vec![150, 200]));
    gene.set_thick_start(Some(119));
    gene.set_thick_end(Some(180));

    let mut buf = Vec::new();
    Writer::<Gtf>::from_record(&gene, &mut buf).unwrap();
    let text = String::from_utf8(buf).unwrap();
    let lines: Vec<&str> = text.trim_end().split('\n').collect();
    assert_eq!(lines.len(), 7);

    assert!(lines[0].starts_with("chr1\tgenepred\ttranscript\t100\t200\t.\t+\t.\t"));
    assert!(lines[1].starts_with("chr1\tgenepred\texon\t100\t150\t.\t+\t.\t"));
    assert!(lines[2].starts_with("chr1\tgenepred\texon\t170\t200\t.\t+\t.\t"));

    let cds_first = lines
        .iter()
        .find(|l| l.contains("\tCDS\t120\t150"))
        .unwrap();
    let cds_second = lines
        .iter()
        .find(|l| l.contains("\tCDS\t170\t180"))
        .unwrap();
    assert!(cds_first.ends_with("\tgene_id \"tx1\"; transcript_id \"tx1\"; tag \"a,b\";"));
    assert!(cds_first.contains("\t+\t0\t"));
    assert!(cds_second.contains("\t+\t2\t"));

    let start_codon = lines.iter().find(|l| l.contains("start_codon")).unwrap();
    assert!(start_codon.contains("\t120\t122\t.\t+\t.\t"));
    let stop_codon = lines.iter().find(|l| l.contains("stop_codon")).unwrap();
    assert!(stop_codon.contains("\t178\t180\t.\t+\t.\t"));
}

#[test]
fn write_gff_reverse_strand_with_phases() {
    let mut gene = GenePred::from_coords(b"chr2".to_vec(), 0, 90, Extras::new());
    gene.set_name(Some(b"txR".to_vec()));
    gene.set_strand(Some(Strand::Reverse));
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![0, 40]));
    gene.set_block_ends(Some(vec![20, 80]));
    gene.set_thick_start(Some(10));
    gene.set_thick_end(Some(80));

    let mut buf = Vec::new();
    Writer::<Gff>::from_record(&gene, &mut buf).unwrap();
    let text = String::from_utf8(buf).unwrap();
    let lines: Vec<&str> = text.trim_end().split('\n').collect();
    assert_eq!(lines.len(), 7);

    let cds_early = lines.iter().find(|l| l.contains("\tCDS\t11\t20")).unwrap();
    assert!(cds_early.contains("\t-\t2\t"));

    let cds_late = lines.iter().find(|l| l.contains("\tCDS\t41\t80")).unwrap();
    assert!(cds_late.contains("\t-\t0\t"));

    let start_codon = lines.iter().find(|l| l.contains("start_codon")).unwrap();
    assert!(start_codon.contains("\t78\t80\t.\t-\t.\t"));

    let stop_codon = lines.iter().find(|l| l.contains("stop_codon")).unwrap();
    assert!(stop_codon.contains("\t11\t13\t.\t-\t.\t"));
}

#[test]
fn write_bed12_preserves_blocks() {
    let mut gene = GenePred::from_coords(b"chr3".to_vec(), 100, 260, Extras::new());
    gene.set_name(Some(b"txBed".to_vec()));
    gene.set_strand(Some(Strand::Forward));
    gene.set_thick_start(Some(120));
    gene.set_thick_end(Some(240));
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![100, 200]));
    gene.set_block_ends(Some(vec![150, 260]));

    let mut buf = Vec::new();
    Writer::<Bed12>::from_record(&gene, &mut buf).unwrap();
    let text = String::from_utf8(buf).unwrap();
    assert_eq!(
        text.trim_end(),
        "chr3\t100\t260\ttxBed\t0\t+\t120\t240\t0,0,0\t2\t50,60,\t0,100,"
    );
}

#[test]
fn write_bed3_orders_numeric_extras() {
    let mut extras = Extras::new();
    extras.insert(b"4".to_vec(), ExtraValue::Scalar(b"fourth".to_vec()));
    extras.insert(b"3".to_vec(), ExtraValue::Scalar(b"third".to_vec()));
    extras.insert(b"note".to_vec(), ExtraValue::Scalar(b"keep".to_vec()));
    let gene = GenePred::from_coords(b"chr4".to_vec(), 10, 20, extras);

    let mut buf = Vec::new();
    let opts = WriterOptions::new().include_non_numeric_extras(true);
    Writer::<Bed3>::from_record_with_options(&gene, &mut buf, &opts).unwrap();
    let text = String::from_utf8(buf).unwrap();
    assert_eq!(text.trim_end(), "chr4\t10\t20\tthird\tfourth\tnote=keep");
}

#[test]
fn write_bed3_skips_non_numeric_by_default() {
    let mut extras = Extras::new();
    extras.insert(b"note".to_vec(), ExtraValue::Scalar(b"keep".to_vec()));
    let gene = GenePred::from_coords(b"chr5".to_vec(), 0, 10, extras);

    let mut buf = Vec::new();
    Writer::<Bed3>::from_record(&gene, &mut buf).unwrap();
    let text = String::from_utf8(buf).unwrap();
    assert_eq!(text.trim_end(), "chr5\t0\t10");
}

#[test]
fn write_gtf_gene_transcript_first() {
    let path = "tests/data/bed12_extra.bed";
    let mut reader: Reader<Bed12> = Reader::from_path_with_additional_fields(path, 2).unwrap();
    let record = reader.records().next().unwrap().unwrap();

    let mut buf = Vec::new();
    Writer::<Gtf>::from_record(&record, &mut buf).unwrap();
    let text = String::from_utf8(buf).unwrap();
    let first_line = text.lines().next().unwrap();
    let attrs = first_line.split('\t').last().unwrap();
    let mut parts = attrs.split(';').filter(|s| !s.trim().is_empty());

    assert!(parts
        .next()
        .unwrap()
        .trim_start()
        .starts_with("gene_id \"txB\""));
    assert!(parts
        .next()
        .unwrap()
        .trim_start()
        .starts_with("transcript_id \"txB\""));
}

#[test]
fn gtf_to_bed_includes_codons_in_cds_bounds() {
    let path = "tests/data/codons.gtf";
    let mut reader: Reader<Gtf> = Reader::from_path(path).unwrap();
    let record = reader.records().next().unwrap().unwrap();

    assert_eq!(record.thick_start(), Some(69));
    assert_eq!(record.thick_end(), Some(200));

    let mut buf = Vec::new();
    Writer::<Bed12>::from_record(&record, &mut buf).unwrap();
    let text = String::from_utf8(buf).unwrap();
    let fields: Vec<&str> = text.trim_end().split('\t').collect();
    assert_eq!(fields[6], "69");
    assert_eq!(fields[7], "200");
}
