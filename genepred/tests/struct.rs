use std::collections::HashMap;

use genepred::bed::{Bed12, Bed3, Bed4, Bed5, Bed6, Bed8, Bed9};
use genepred::{
    BedFormat, ExtraValue, Extras, GeneLine, GenePred, Gff, Gtf, GxfOptions, Reader, ReaderOptions,
    Strand, TranscriptParent, Writer, WriterOptions,
};

#[test]
fn test_genepred_from_coords() {
    let mut extras = Extras::new();
    extras.insert(b"key1".to_vec(), ExtraValue::Scalar(b"extra1".to_vec()));
    extras.insert(b"key2".to_vec(), ExtraValue::Scalar(b"extra2".to_vec()));
    let gene = GenePred::from_coords(b"chr1".to_vec(), 10, 20, extras.clone());
    assert_eq!(gene.chrom(), b"chr1".as_ref());
    assert_eq!(gene.start(), 10);
    assert_eq!(gene.end(), 20);
    assert_eq!(gene.extras(), &extras);
    assert!(gene.name().is_none());
}

#[test]
fn test_genepred_from_bed3() {
    let bed3 = Bed3 {
        chrom: b"chr1".to_vec(),
        start: 10,
        end: 20,
        extras: Extras::new(),
    };
    let gene: GenePred = bed3.into();
    assert_eq!(gene.chrom(), b"chr1".as_ref());
    assert_eq!(gene.start(), 10);
    assert_eq!(gene.end(), 20);
    assert!(gene.name().is_none());
}

#[test]
fn test_genepred_from_bed4() {
    let bed4 = Bed4 {
        chrom: b"chr1".to_vec(),
        start: 10,
        end: 20,
        name: b"geneA".to_vec(),
        extras: Extras::new(),
    };
    let gene: GenePred = bed4.into();
    assert_eq!(gene.chrom(), b"chr1".as_ref());
    assert_eq!(gene.start(), 10);
    assert_eq!(gene.end(), 20);
    assert_eq!(gene.name().unwrap(), b"geneA".as_ref());
}

#[test]
fn test_genepred_from_bed6() {
    let bed6 = Bed6 {
        chrom: b"chr1".to_vec(),
        start: 10,
        end: 20,
        name: b"geneA".to_vec(),
        score: 100,
        strand: Strand::Forward,
        extras: Extras::new(),
    };
    let gene: GenePred = bed6.into();
    assert_eq!(gene.chrom(), b"chr1".as_ref());
    assert_eq!(gene.start(), 10);
    assert_eq!(gene.end(), 20);
    assert_eq!(gene.name().unwrap(), b"geneA".as_ref());
    assert_eq!(gene.strand().unwrap(), Strand::Forward);
}

#[test]
fn test_genepred_from_bed12() {
    let bed12 = Bed12 {
        chrom: b"chr1".to_vec(),
        start: 10,
        end: 100,
        name: b"geneA".to_vec(),
        score: 1000,
        strand: Strand::Forward,
        thick_start: 10,
        thick_end: 100,
        item_rgb: genepred::bed::Rgb(255, 0, 0),
        block_count: 2,
        block_sizes: vec![10, 20],
        block_starts: vec![0, 30],
        extras: Extras::new(),
    };
    let gene: GenePred = bed12.into();
    assert_eq!(gene.chrom(), b"chr1".as_ref());
    assert_eq!(gene.start(), 10);
    assert_eq!(gene.end(), 100);
    assert_eq!(gene.name().unwrap(), b"geneA".as_ref());
    assert_eq!(gene.strand().unwrap(), Strand::Forward);
    assert_eq!(gene.thick_start().unwrap(), 10);
    assert_eq!(gene.thick_end().unwrap(), 100);
    assert_eq!(gene.block_count().unwrap(), 2);
    assert_eq!(gene.block_starts().unwrap(), &[10u64, 40]);
    assert_eq!(gene.block_ends().unwrap(), &[20u64, 60]);
}

#[test]
fn test_genepred_len_is_empty() {
    let gene = GenePred::from_coords(b"chr1".to_vec(), 10, 20, Extras::new());
    assert_eq!(gene.len(), 10);
    assert!(!gene.is_empty());

    let empty_gene = GenePred::from_coords(b"chr1".to_vec(), 10, 10, Extras::new());
    assert_eq!(empty_gene.len(), 0);
    assert!(empty_gene.is_empty());

    let inverted_gene = GenePred::from_coords(b"chr1".to_vec(), 20, 10, Extras::new());
    assert_eq!(inverted_gene.len(), 0);
    assert!(inverted_gene.is_empty());
}

#[test]
fn test_genepred_setters() {
    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 20, Extras::new());
    gene.set_chrom(b"chrX".to_vec());
    gene.set_start(5);
    gene.set_end(25);
    gene.set_name(Some(b"new_gene".to_vec()));
    gene.set_strand(Some(Strand::Reverse));
    gene.set_thick_start(Some(6));
    gene.set_thick_end(Some(24));
    gene.set_block_count(Some(1));
    gene.set_block_starts(Some(vec![5]));
    gene.set_block_ends(Some(vec![25]));
    let mut extras_map = Extras::new();
    let extra_key = b"extra".to_vec();
    extras_map.insert(extra_key.clone(), ExtraValue::Scalar(b"new_extra".to_vec()));
    gene.set_extras(extras_map.clone());

    assert_eq!(gene.chrom(), b"chrX".as_ref());
    assert_eq!(gene.start(), 5);
    assert_eq!(gene.end(), 25);
    assert_eq!(gene.name().unwrap(), b"new_gene".as_ref());
    assert_eq!(gene.strand().unwrap(), Strand::Reverse);
    assert_eq!(gene.thick_start().unwrap(), 6);
    assert_eq!(gene.thick_end().unwrap(), 24);
    assert_eq!(gene.block_count().unwrap(), 1);
    assert_eq!(gene.block_starts().unwrap(), &[5]);
    assert_eq!(gene.block_ends().unwrap(), &[25]);
    assert_eq!(gene.extras(), &extras_map);

    gene.add_extra(extra_key.clone(), b"another_extra".to_vec());
    if let Some(value) = extras_map.get_mut(&extra_key) {
        value.push(b"another_extra".to_vec());
    }
    assert_eq!(gene.extras(), &extras_map);
    gene.clear_extras();
    assert!(gene.extras().is_empty());
}

#[test]
fn test_genepred_exons() {
    // No blocks defined
    let gene1 = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    assert_eq!(gene1.exons(), vec![(10, 100)]);

    // With blocks defined
    let mut gene2 = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    gene2.set_block_count(Some(2));
    gene2.set_block_starts(Some(vec![10, 40]));
    gene2.set_block_ends(Some(vec![20, 60])); // Absolute coordinates
    assert_eq!(gene2.exons(), vec![(10, 20), (40, 60)]);

    // With empty blocks
    let mut gene3 = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    gene3.set_block_count(Some(0));
    gene3.set_block_starts(Some(vec![]));
    gene3.set_block_ends(Some(vec![]));
    assert_eq!(gene3.exons(), vec![(10, 100)]);
}

#[test]
fn test_genepred_introns() {
    // No introns (single exon)
    let gene1 = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    assert!(gene1.introns().is_empty());

    // With multiple exons
    let mut gene2 = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    gene2.set_block_count(Some(2));
    gene2.set_block_starts(Some(vec![10, 40])); // Exons: (10,20), (40,60)
    gene2.set_block_ends(Some(vec![20, 60]));
    assert_eq!(gene2.introns(), vec![(20, 40)]);

    // With more exons
    let mut gene3 = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    gene3.set_block_count(Some(3));
    gene3.set_block_starts(Some(vec![10, 30, 50])); // Exons: (10,20), (30,40), (50,60)
    gene3.set_block_ends(Some(vec![20, 40, 60]));
    assert_eq!(gene3.introns(), vec![(20, 30), (40, 50)]);
}

#[test]
fn test_genepred_exonic_intronic_length() {
    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![10, 40])); // Exons: (10,20), (40,60)
    gene.set_block_ends(Some(vec![20, 60]));
    assert_eq!(gene.exonic_length(), 30);
    assert_eq!(gene.intronic_length(), 20);

    let gene_no_blocks = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    assert_eq!(gene_no_blocks.exonic_length(), 90); // (100-10)
    assert_eq!(gene_no_blocks.intronic_length(), 0);
}

#[test]
fn test_genepred_coding_exons_cds_length() {
    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![10, 40])); // Exons: (10,20), (40,60)
    gene.set_block_ends(Some(vec![20, 60]));

    // No thick regions
    assert!(gene.coding_exons().is_empty());
    assert_eq!(gene.cds_length(), 0);

    // Thick regions overlapping exons
    gene.set_thick_start(Some(15));
    gene.set_thick_end(Some(50)); // Overlaps (10,20) -> (15,20), Overlaps (40,60) -> (40,50)
    assert_eq!(gene.coding_exons(), vec![(15, 20), (40, 50)]);
    assert_eq!(gene.cds_length(), 5 + 10); // 15

    // Thick regions not overlapping exons
    gene.set_thick_start(Some(70));
    gene.set_thick_end(Some(80));
    assert!(gene.coding_exons().is_empty());
    assert_eq!(gene.cds_length(), 0);
}

#[test]
fn test_extra_value_conversion_and_empty_helpers() {
    let scalar = ExtraValue::Scalar(b"value1".to_vec());
    assert_eq!(scalar.clone().into_inner(), vec![b"value1".to_vec()]);
    assert_eq!(scalar.clone().into_scalar(), Some(b"value1".to_vec()));
    assert_eq!(scalar.into_array(), None);

    let array = ExtraValue::Array(vec![b"value1".to_vec(), b"value2".to_vec()]);
    assert_eq!(
        array.clone().into_inner(),
        vec![b"value1".to_vec(), b"value2".to_vec()]
    );
    assert_eq!(array.clone().into_scalar(), None);
    assert_eq!(
        array.into_array(),
        Some(vec![b"value1".to_vec(), b"value2".to_vec()])
    );

    assert!(ExtraValue::Scalar(Vec::new()).is_empty());
    assert!(ExtraValue::Array(Vec::new()).is_empty());
    assert!(!ExtraValue::Scalar(b"v".to_vec()).is_empty());
    assert!(!ExtraValue::Array(vec![Vec::new()]).is_empty());
}

#[test]
fn test_genepred_get_extra() {
    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 20, Extras::new());
    gene.add_extra("gene_id", "gene1");
    gene.add_extra("gene_id", "gene2");

    assert_eq!(gene.get_extra(b"missing"), None);
    assert_eq!(
        gene.get_extra(b"gene_id"),
        Some(&ExtraValue::Array(vec![
            b"gene1".to_vec(),
            b"gene2".to_vec()
        ]))
    );
}

#[test]
fn test_genepred_utr_exons_and_utr_length() {
    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 120, Extras::new());
    gene.set_block_count(Some(3));
    gene.set_block_starts(Some(vec![10, 40, 100]));
    gene.set_block_ends(Some(vec![30, 90, 120]));
    gene.set_thick_start(Some(50));
    gene.set_thick_end(Some(80));

    assert_eq!(
        gene.utr_exons(),
        vec![(10, 30), (40, 50), (80, 90), (100, 120)]
    );
    assert_eq!(gene.utr_length(), 60);
}

#[test]
fn test_genepred_utr_exons_include_boundary_touching_exons() {
    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 40, Extras::new());
    gene.set_block_count(Some(3));
    gene.set_block_starts(Some(vec![10, 20, 30]));
    gene.set_block_ends(Some(vec![20, 30, 40]));
    gene.set_thick_start(Some(20));
    gene.set_thick_end(Some(30));

    assert_eq!(gene.utr_exons(), vec![(10, 20), (30, 40)]);
}

#[test]
fn test_genepred_strand_aware_utrs() {
    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 120, Extras::new());
    gene.set_block_count(Some(3));
    gene.set_block_starts(Some(vec![10, 40, 100]));
    gene.set_block_ends(Some(vec![30, 90, 120]));
    gene.set_thick_start(Some(50));
    gene.set_thick_end(Some(80));

    gene.set_strand(Some(Strand::Forward));
    assert_eq!(gene.five_prime_utr(), vec![(10, 30), (40, 50)]);
    assert_eq!(gene.three_prime_utr(), vec![(80, 90), (100, 120)]);

    gene.set_strand(Some(Strand::Reverse));
    assert_eq!(gene.five_prime_utr(), vec![(80, 90), (100, 120)]);
    assert_eq!(gene.three_prime_utr(), vec![(10, 30), (40, 50)]);

    gene.set_strand(Some(Strand::Unknown));
    assert!(gene.five_prime_utr().is_empty());
    assert!(gene.three_prime_utr().is_empty());
}

#[test]
fn test_genepred_unnest_extras() {
    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 20, Extras::new());
    gene.add_extra(b"attrs".to_vec(), b"tag1=value1;tag2=value2".to_vec());
    gene.add_extra(b"attrs".to_vec(), b"tag3=value3".to_vec());

    assert_eq!(
        gene.unnest_extras(";"),
        vec![
            b"tag1=value1".to_vec(),
            b"tag2=value2".to_vec(),
            b"tag3=value3".to_vec()
        ]
    );

    assert_eq!(
        gene.unnest_extras("="),
        vec![
            b"tag1".to_vec(),
            b"value1;tag2".to_vec(),
            b"value2".to_vec(),
            b"tag3".to_vec(),
            b"value3".to_vec()
        ]
    );
}

#[test]
fn test_genepred_overlaps() {
    let gene = GenePred::from_coords(b"chr1".to_vec(), 10, 20, Extras::new());

    // Complete overlap
    assert!(gene.overlaps(5, 25));
    // Partial overlap (start)
    assert!(gene.overlaps(5, 15));
    // Partial overlap (end)
    assert!(gene.overlaps(15, 25));
    // Exact overlap
    assert!(gene.overlaps(10, 20));
    // No overlap (before)
    assert!(!gene.overlaps(0, 5));
    // No overlap (after)
    assert!(!gene.overlaps(25, 30));
    // Touch (start)
    assert!(!gene.overlaps(5, 10));
    // Touch (end)
    assert!(!gene.overlaps(20, 25));
}

#[test]
fn test_genepred_exon_overlaps() {
    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![10, 40])); // Exons: (10,20), (40,60)
    gene.set_block_ends(Some(vec![20, 60]));

    // Overlaps first exon
    assert!(gene.exon_overlaps(5, 15));
    // Overlaps second exon
    assert!(gene.exon_overlaps(45, 65));
    // Overlaps both exons (and intron)
    assert!(gene.exon_overlaps(15, 50));
    // No overlap
    assert!(!gene.exon_overlaps(25, 35));
    assert!(!gene.exon_overlaps(0, 5));
    assert!(!gene.exon_overlaps(70, 80));
}

#[test]
fn test_genepred_exon_intron_count() {
    let gene1 = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    assert_eq!(gene1.exon_count(), 1);
    assert_eq!(gene1.intron_count(), 0);

    let mut gene2 = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    gene2.set_block_count(Some(2));
    gene2.set_block_starts(Some(vec![10, 40])); // Exons: (10,20), (40,60)
    gene2.set_block_ends(Some(vec![20, 60]));
    assert_eq!(gene2.exon_count(), 2);
    assert_eq!(gene2.intron_count(), 1);

    let mut gene3 = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    gene3.set_block_count(Some(0));
    gene3.set_block_starts(Some(vec![]));
    gene3.set_block_ends(Some(vec![]));
    assert_eq!(gene3.exon_count(), 1);
    assert_eq!(gene3.intron_count(), 0);
}

#[test]
fn test_genepred_display() {
    let gene1 = GenePred::from_coords(b"chr1".to_vec(), 10, 20, Extras::new());
    assert_eq!(format!("{}", gene1), ".\tchr1\t.\t10\t20\t10\t20\t0\t\t");

    let mut gene2 = GenePred::from_coords(b"chr1".to_vec(), 10, 20, Extras::new());
    gene2.set_name(Some(b"geneA".to_vec()));
    gene2.set_strand(Some(Strand::Forward));
    assert_eq!(
        format!("{}", gene2),
        "geneA\tchr1\t+\t10\t20\t10\t20\t0\t\t"
    );

    let mut gene3 = GenePred::from_coords(b"chr1".to_vec(), 10, 100, Extras::new());
    gene3.set_name(Some(b"geneB".to_vec()));
    gene3.set_strand(Some(Strand::Reverse));
    gene3.set_thick_start(Some(10));
    gene3.set_thick_end(Some(100));
    gene3.set_block_count(Some(2));
    gene3.set_block_starts(Some(vec![10, 40]));
    gene3.set_block_ends(Some(vec![20, 60]));
    let mut extras = Extras::new();
    extras.insert(b"first".to_vec(), ExtraValue::Scalar(b"extra1".to_vec()));
    extras.insert(b"second".to_vec(), ExtraValue::Scalar(b"extra2".to_vec()));
    gene3.set_extras(extras);
    assert_eq!(
        format!("{}", gene3),
        "geneB\tchr1\t-\t10\t100\t10\t100\t2\t10,40\t20,60\tfirst=extra1\tsecond=extra2"
    );
}

#[test]
fn test_genepred_to_bed_layouts() {
    let mut extras = Extras::new();
    extras.insert(b"4".to_vec(), ExtraValue::Scalar(b"ignored".to_vec()));
    extras.insert(b"13".to_vec(), ExtraValue::Scalar(b"ignored_too".to_vec()));

    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 100, extras);
    gene.set_name(Some(b"geneA".to_vec()));
    gene.set_strand(Some(Strand::Forward));
    gene.set_thick_start(Some(20));
    gene.set_thick_end(Some(90));
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![10, 40]));
    gene.set_block_ends(Some(vec![20, 60]));

    assert_eq!(
        String::from_utf8(gene.to_bed::<Bed3>()).unwrap(),
        "chr1\t10\t100"
    );
    assert_eq!(
        String::from_utf8(gene.to_bed::<Bed4>()).unwrap(),
        "chr1\t10\t100\tgeneA"
    );
    assert_eq!(
        String::from_utf8(gene.to_bed::<Bed5>()).unwrap(),
        "chr1\t10\t100\tgeneA\t0"
    );
    assert_eq!(
        String::from_utf8(gene.to_bed::<Bed6>()).unwrap(),
        "chr1\t10\t100\tgeneA\t0\t+"
    );
    assert_eq!(
        String::from_utf8(gene.to_bed::<Bed8>()).unwrap(),
        "chr1\t10\t100\tgeneA\t0\t+\t20\t90"
    );
    assert_eq!(
        String::from_utf8(gene.to_bed::<Bed9>()).unwrap(),
        "chr1\t10\t100\tgeneA\t0\t+\t20\t90\t0,0,0"
    );
    assert_eq!(
        String::from_utf8(gene.to_bed::<Bed12>()).unwrap(),
        "chr1\t10\t100\tgeneA\t0\t+\t20\t90\t0,0,0\t2\t10,20,\t0,30,"
    );
}

#[test]
fn test_genepred_to_bed_with_additional_fields() {
    let mut extras = Extras::new();
    extras.insert(b"13".to_vec(), ExtraValue::Scalar(b"foo".to_vec()));
    extras.insert(
        b"14".to_vec(),
        ExtraValue::Array(vec![b"bar".to_vec(), b"baz".to_vec()]),
    );

    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 100, extras);
    gene.set_name(Some(b"geneA".to_vec()));
    gene.set_strand(Some(Strand::Forward));
    gene.set_thick_start(Some(20));
    gene.set_thick_end(Some(90));
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![10, 40]));
    gene.set_block_ends(Some(vec![20, 60]));

    assert_eq!(
        String::from_utf8(gene.to_bed_with_additional_fields::<Bed12>(2)).unwrap(),
        "chr1\t10\t100\tgeneA\t0\t+\t20\t90\t0,0,0\t2\t10,20,\t0,30,\tfoo\tbar,baz"
    );
}

#[test]
#[should_panic(expected = "missing additional BED field key '14'")]
fn test_genepred_to_bed_with_additional_fields_panics_when_missing() {
    let mut extras = Extras::new();
    extras.insert(b"13".to_vec(), ExtraValue::Scalar(b"foo".to_vec()));

    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 100, extras);
    gene.set_name(Some(b"geneA".to_vec()));
    gene.set_strand(Some(Strand::Forward));
    gene.set_thick_start(Some(20));
    gene.set_thick_end(Some(90));
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![10, 40]));
    gene.set_block_ends(Some(vec![20, 60]));

    let _ = gene.to_bed_with_additional_fields::<Bed12>(2);
}

#[test]
fn test_genepred_to_gtf_with_mapping_and_additional_fields() {
    let mut extras = Extras::new();
    extras.insert(b"13".to_vec(), ExtraValue::Scalar(b"foo".to_vec()));
    extras.insert(
        b"14".to_vec(),
        ExtraValue::Array(vec![b"bar".to_vec(), b"baz".to_vec()]),
    );

    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 99, 200, extras);
    gene.set_name(Some(b"tx1".to_vec()));
    gene.set_strand(Some(Strand::Forward));
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![99, 169]));
    gene.set_block_ends(Some(vec![150, 200]));
    gene.set_thick_start(Some(119));
    gene.set_thick_end(Some(180));

    let mut mapping = HashMap::new();
    mapping.insert("tx1".to_string(), "gene1".to_string());

    let lines = gene.to_gxf_with_additional_fields::<Gtf>(2, Some(&mapping));
    let text: Vec<String> = lines
        .into_iter()
        .map(|line| String::from_utf8(line).unwrap())
        .collect();

    assert_eq!(text.len(), 8);
    assert_eq!(
        text[0],
        "chr1\tgenepred\tgene\t100\t200\t.\t+\t.\tgene_id \"gene1\"; 13 \"foo\"; 14 \"bar,baz\";"
    );
    assert_eq!(
        text[1],
        "chr1\tgenepred\ttranscript\t100\t200\t.\t+\t.\tgene_id \"gene1\"; transcript_id \"tx1\"; 13 \"foo\"; 14 \"bar,baz\";"
    );
    assert_eq!(
        text[2],
        "chr1\tgenepred\texon\t100\t150\t.\t+\t.\tgene_id \"gene1\"; transcript_id \"tx1\"; exon_number \"1\"; 13 \"foo\"; 14 \"bar,baz\";"
    );
    assert_eq!(
        text[3],
        "chr1\tgenepred\texon\t170\t200\t.\t+\t.\tgene_id \"gene1\"; transcript_id \"tx1\"; exon_number \"2\"; 13 \"foo\"; 14 \"bar,baz\";"
    );
    assert!(text.contains(
        &"chr1\tgenepred\tCDS\t120\t150\t.\t+\t0\tgene_id \"gene1\"; transcript_id \"tx1\"; exon_number \"1\"; 13 \"foo\"; 14 \"bar,baz\";".to_string()
    ));
    assert!(text.contains(
        &"chr1\tgenepred\tCDS\t170\t180\t.\t+\t2\tgene_id \"gene1\"; transcript_id \"tx1\"; exon_number \"2\"; 13 \"foo\"; 14 \"bar,baz\";".to_string()
    ));
    assert!(text.contains(
        &"chr1\tgenepred\tstart_codon\t120\t122\t.\t+\t.\tgene_id \"gene1\"; transcript_id \"tx1\"; exon_number \"1\"; 13 \"foo\"; 14 \"bar,baz\";".to_string()
    ));
    assert!(text.contains(
        &"chr1\tgenepred\tstop_codon\t178\t180\t.\t+\t.\tgene_id \"gene1\"; transcript_id \"tx1\"; exon_number \"2\"; 13 \"foo\"; 14 \"bar,baz\";".to_string()
    ));
}

#[test]
fn test_genepred_to_gff_standard_hierarchy() {
    let mut extras = Extras::new();
    extras.insert(b"13".to_vec(), ExtraValue::Scalar(b"foo".to_vec()));

    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 99, 200, extras);
    gene.set_name(Some(b"tx1".to_vec()));
    gene.set_strand(Some(Strand::Forward));
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![99, 169]));
    gene.set_block_ends(Some(vec![150, 200]));
    gene.set_thick_start(Some(119));
    gene.set_thick_end(Some(180));

    let mut mapping = HashMap::new();
    mapping.insert("tx1".to_string(), "gene1".to_string());

    let lines = gene.to_gxf_with_additional_fields::<Gff>(1, Some(&mapping));
    let text: Vec<String> = lines
        .into_iter()
        .map(|line| String::from_utf8(line).unwrap())
        .collect();

    assert_eq!(
        text[0],
        "chr1\tgenepred\tgene\t100\t200\t.\t+\t.\tID=gene1;13=foo;"
    );
    assert_eq!(
        text[1],
        "chr1\tgenepred\tmRNA\t100\t200\t.\t+\t.\tID=tx1;Parent=gene1;13=foo;"
    );
    assert_eq!(
        text[2],
        "chr1\tgenepred\texon\t100\t150\t.\t+\t.\tParent=tx1;exon_number=1;13=foo;"
    );
    assert_eq!(
        text[3],
        "chr1\tgenepred\texon\t170\t200\t.\t+\t.\tParent=tx1;exon_number=2;13=foo;"
    );
    assert!(text.contains(
        &"chr1\tgenepred\tCDS\t120\t150\t.\t+\t0\tParent=tx1;exon_number=1;13=foo;".to_string()
    ));
    assert!(text.contains(
        &"chr1\tgenepred\tCDS\t170\t180\t.\t+\t2\tParent=tx1;exon_number=2;13=foo;".to_string()
    ));
    assert!(text.contains(
        &"chr1\tgenepred\tstart_codon\t120\t122\t.\t+\t.\tParent=tx1;exon_number=1;13=foo;"
            .to_string()
    ));
    assert!(text.contains(
        &"chr1\tgenepred\tstop_codon\t178\t180\t.\t+\t.\tParent=tx1;exon_number=2;13=foo;"
            .to_string()
    ));
}

#[test]
fn test_genepred_to_gxf_mapping_uses_name_key() {
    let mut extras = Extras::new();
    extras.insert(
        b"transcript_id".to_vec(),
        ExtraValue::Scalar(b"txCanonical".to_vec()),
    );

    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 10, 20, extras);
    gene.set_name(Some(b"txAlias".to_vec()));

    let mut mapping = HashMap::new();
    mapping.insert("txAlias".to_string(), "geneAlias".to_string());

    let lines = gene.to_gxf::<Gtf>(Some(&mapping));
    let transcript = String::from_utf8(lines[1].clone()).unwrap();

    assert_eq!(
        transcript,
        "chr1\tgenepred\ttranscript\t11\t20\t.\t.\t.\tgene_id \"geneAlias\"; transcript_id \"txCanonical\";"
    );
}

#[test]
fn test_genepred_to_gxf_uses_strand_aware_exon_numbers() {
    let mut gene = GenePred::from_coords(b"chr2".to_vec(), 0, 90, Extras::new());
    gene.set_name(Some(b"txR".to_vec()));
    gene.set_strand(Some(Strand::Reverse));
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![0, 40]));
    gene.set_block_ends(Some(vec![20, 80]));
    gene.set_thick_start(Some(10));
    gene.set_thick_end(Some(80));

    let lines = gene.to_gxf::<Gff>(None);
    let text: Vec<String> = lines
        .into_iter()
        .map(|line| String::from_utf8(line).unwrap())
        .collect();

    assert_eq!(
        text[2],
        "chr2\tgenepred\texon\t1\t20\t.\t-\t.\tParent=txR;exon_number=2;"
    );
    assert_eq!(
        text[3],
        "chr2\tgenepred\texon\t41\t80\t.\t-\t.\tParent=txR;exon_number=1;"
    );
}

#[test]
#[should_panic(expected = "unsupported GXF layout")]
fn test_genepred_to_gxf_panics_for_bed_layout() {
    let gene = GenePred::from_coords(b"chr1".to_vec(), 10, 20, Extras::new());
    let _ = gene.to_gxf::<Bed12>(None);
}

#[test]
#[should_panic(expected = "missing numeric GXF additional fields")]
fn test_genepred_to_gxf_with_additional_fields_panics_when_missing_numeric_extras() {
    let mut extras = Extras::new();
    extras.insert(b"13".to_vec(), ExtraValue::Scalar(b"foo".to_vec()));
    extras.insert(b"note".to_vec(), ExtraValue::Scalar(b"keep".to_vec()));

    let gene = GenePred::from_coords(b"chr1".to_vec(), 10, 20, extras);
    let _ = gene.to_gxf_with_additional_fields::<Gtf>(2, None);
}

/// A multi-exon, coding `+` strand gene with two numeric extras (`13`, `14`),
/// matching the fixtures used by the standard GTF/GFF hierarchy tests.
fn standard_coding_gene() -> GenePred {
    let mut extras = Extras::new();
    extras.insert(b"13".to_vec(), ExtraValue::Scalar(b"foo".to_vec()));
    extras.insert(
        b"14".to_vec(),
        ExtraValue::Array(vec![b"bar".to_vec(), b"baz".to_vec()]),
    );

    let mut gene = GenePred::from_coords(b"chr1".to_vec(), 99, 200, extras);
    gene.set_name(Some(b"tx1".to_vec()));
    gene.set_strand(Some(Strand::Forward));
    gene.set_block_count(Some(2));
    gene.set_block_starts(Some(vec![99, 169]));
    gene.set_block_ends(Some(vec![150, 200]));
    gene.set_thick_start(Some(119));
    gene.set_thick_end(Some(180));
    gene
}

fn mapping_tx1_gene1() -> HashMap<String, String> {
    let mut mapping = HashMap::new();
    mapping.insert("tx1".to_string(), "gene1".to_string());
    mapping
}

fn as_text(lines: Vec<Vec<u8>>) -> Vec<String> {
    lines
        .into_iter()
        .map(|line| String::from_utf8(line).unwrap())
        .collect()
}

#[test]
fn test_to_gxf_with_options_include_matches_legacy_methods() {
    let gene = standard_coding_gene();
    let mapping = mapping_tx1_gene1();

    assert_eq!(
        gene.to_gxf_with_options::<Gtf>(&GxfOptions::default(), Some(&mapping)),
        gene.to_gxf::<Gtf>(Some(&mapping))
    );

    let opts = GxfOptions {
        additional_fields: 2,
        gene_line: GeneLine::Include,
        ..Default::default()
    };
    assert_eq!(
        gene.to_gxf_with_options::<Gtf>(&opts, Some(&mapping)),
        gene.to_gxf_with_additional_fields::<Gtf>(2, Some(&mapping))
    );
}

#[test]
fn test_to_gxf_omit_drops_only_gene_line_gtf() {
    let gene = standard_coding_gene();
    let mapping = mapping_tx1_gene1();
    let include = GxfOptions {
        additional_fields: 2,
        gene_line: GeneLine::Include,
        ..Default::default()
    };
    let omit = GxfOptions {
        additional_fields: 2,
        gene_line: GeneLine::Omit,
        ..Default::default()
    };

    let included = gene.to_gxf_with_options::<Gtf>(&include, Some(&mapping));
    let omitted = gene.to_gxf_with_options::<Gtf>(&omit, Some(&mapping));

    // Omit is exactly Include without its leading `gene` line.
    assert_eq!(omitted.len() + 1, included.len());
    assert_eq!(omitted, included[1..].to_vec());

    let text = as_text(omitted);
    assert!(text[0].contains("\ttranscript\t"));
    assert!(text
        .iter()
        .all(|line| line.split('\t').nth(2) != Some("gene")));
}

#[test]
fn test_to_gxf_omit_drops_only_gene_line_gff() {
    let gene = standard_coding_gene();
    let mapping = mapping_tx1_gene1();
    let include = GxfOptions {
        additional_fields: 1,
        gene_line: GeneLine::Include,
        ..Default::default()
    };
    let omit = GxfOptions {
        additional_fields: 1,
        gene_line: GeneLine::Omit,
        ..Default::default()
    };

    let included = gene.to_gxf_with_options::<Gff>(&include, Some(&mapping));
    let omitted = gene.to_gxf_with_options::<Gff>(&omit, Some(&mapping));

    assert_eq!(omitted, included[1..].to_vec());

    let text = as_text(omitted);
    assert!(text[0].contains("\tmRNA\t"));
    assert!(text
        .iter()
        .all(|line| line.split('\t').nth(2) != Some("gene")));
}

#[test]
fn test_to_gxf_gene_line_gtf() {
    let gene = standard_coding_gene();
    let mapping = mapping_tx1_gene1();

    let line = String::from_utf8(gene.to_gxf_gene_line::<Gtf>(2, Some(&mapping))).unwrap();
    assert_eq!(
        line,
        "chr1\tgenepred\tgene\t100\t200\t.\t+\t.\tgene_id \"gene1\"; 13 \"foo\"; 14 \"bar,baz\";"
    );
}

#[test]
fn test_to_gxf_gene_line_gff() {
    let gene = standard_coding_gene();
    let mapping = mapping_tx1_gene1();

    let line = String::from_utf8(gene.to_gxf_gene_line::<Gff>(1, Some(&mapping))).unwrap();
    assert_eq!(
        line,
        "chr1\tgenepred\tgene\t100\t200\t.\t+\t.\tID=gene1;13=foo;"
    );
}

#[test]
fn test_to_gxf_gene_line_reflects_record_span() {
    // Aggregating callers pass a record carrying the union-of-isoforms span.
    let mut gene = GenePred::from_coords(b"chrZ".to_vec(), 900, 2200, Extras::new());
    gene.set_name(Some(b"txU".to_vec()));
    gene.set_strand(Some(Strand::Reverse));

    let line = String::from_utf8(gene.to_gxf_gene_line::<Gtf>(0, None)).unwrap();
    assert_eq!(
        line,
        "chrZ\tgenepred\tgene\t901\t2200\t.\t-\t.\tgene_id \"txU\";"
    );
}

#[test]
#[should_panic(expected = "unsupported GXF layout")]
fn test_to_gxf_gene_line_panics_for_bed_layout() {
    let gene = GenePred::from_coords(b"chr1".to_vec(), 10, 20, Extras::new());
    let _ = gene.to_gxf_gene_line::<Bed12>(0, None);
}

/// Without a gene mapping the gene id falls back to the transcript id. In GFF
/// that would make the `gene` and `mRNA` rows collide on `ID`, so the gene id
/// is prefixed with `gene-`; the transcript `Parent` then points at the gene,
/// never at itself.
#[test]
fn test_to_gxf_no_mapping_gff_synthesizes_distinct_gene_id() {
    let gene = standard_coding_gene();
    let text = as_text(gene.to_gxf::<Gff>(None));

    assert_eq!(
        text[0],
        "chr1\tgenepred\tgene\t100\t200\t.\t+\t.\tID=gene-tx1;"
    );
    assert_eq!(
        text[1],
        "chr1\tgenepred\tmRNA\t100\t200\t.\t+\t.\tID=tx1;Parent=gene-tx1;"
    );
    // Children still hang off the transcript, not the synthesized gene.
    assert_eq!(
        text[2],
        "chr1\tgenepred\texon\t100\t150\t.\t+\t.\tParent=tx1;exon_number=1;"
    );

    // No two features share an `ID` (GFF3 uniqueness), and the transcript is
    // not its own parent.
    let ids: Vec<&str> = text
        .iter()
        .filter_map(|line| {
            line.split('\t')
                .nth(8)?
                .split(';')
                .find_map(|attr| attr.strip_prefix("ID="))
        })
        .collect();
    assert_eq!(ids, vec!["gene-tx1", "tx1"]);
}

/// GTF has no `ID` uniqueness rule, so the bare `gene_id == transcript_id`
/// fallback is left untouched: the `gene-` synthesis is GFF-only.
#[test]
fn test_to_gxf_no_mapping_gtf_keeps_bare_gene_id() {
    let gene = standard_coding_gene();
    let text = as_text(gene.to_gxf::<Gtf>(None));

    assert_eq!(
        text[0],
        "chr1\tgenepred\tgene\t100\t200\t.\t+\t.\tgene_id \"tx1\";"
    );
    assert_eq!(
        text[1],
        "chr1\tgenepred\ttranscript\t100\t200\t.\t+\t.\tgene_id \"tx1\"; transcript_id \"tx1\";"
    );
}

/// No-gene mode (`GeneLine::Omit` + `TranscriptParent::Omit`) makes the
/// transcript top-level: no gene line and no `Parent`, regardless of whether a
/// mapping exists. This is the `--no-gene` contract.
#[test]
fn test_to_gxf_no_gene_mode_makes_transcript_top_level_gff() {
    let gene = standard_coding_gene();
    let mapping = mapping_tx1_gene1();
    let opts = GxfOptions {
        gene_line: GeneLine::Omit,
        transcript_parent: TranscriptParent::Omit,
        ..Default::default()
    };

    let text = as_text(gene.to_gxf_with_options::<Gff>(&opts, Some(&mapping)));

    // First row is the transcript (no gene line) and carries no `Parent`,
    // even though tx1 maps to gene1.
    assert_eq!(text[0], "chr1\tgenepred\tmRNA\t100\t200\t.\t+\t.\tID=tx1;");
    assert!(text
        .iter()
        .all(|line| line.split('\t').nth(2) != Some("gene")));
    assert!(
        !text[0].contains("Parent="),
        "top-level transcript must not carry a Parent: {}",
        text[0]
    );
    // Children still reference the transcript.
    assert_eq!(
        text[1],
        "chr1\tgenepred\texon\t100\t150\t.\t+\t.\tParent=tx1;exon_number=1;"
    );
}

/// Regression for the original `--no-gene` report: with no mapping the gene id
/// equals the transcript id, and emitting `Parent=<gene_id>` produced
/// `Parent=tx1` on `ID=tx1` — a self-referential mRNA that loops GFF parsers.
/// In no-gene mode the transcript is top-level and must have no `Parent`.
#[test]
fn test_to_gxf_no_gene_mode_no_mapping_has_no_self_parent_gff() {
    let gene = standard_coding_gene();
    let opts = GxfOptions {
        gene_line: GeneLine::Omit,
        transcript_parent: TranscriptParent::Omit,
        ..Default::default()
    };

    let text = as_text(gene.to_gxf_with_options::<Gff>(&opts, None));

    assert_eq!(text[0], "chr1\tgenepred\tmRNA\t100\t200\t.\t+\t.\tID=tx1;");
    assert!(
        !text[0].contains("Parent=tx1"),
        "transcript must not be its own parent: {}",
        text[0]
    );
}

/// The release's writer options (release_design.md sec.3), replayed on the +/-
/// mirrored pair of `tests/data/bed12_mirror.bed`. Expected rows are derived by
/// hand below; every row of one record is the mirror of the other's.
///
/// Derivation (0-based, half-open; the records share the span [100, 380) and the
/// mirror map is x -> 480 - x):
///
/// - txMirP (+) exons [100,130) [180,220) [320,380), thick [140,361) -> coding
///   pieces [180,220) 40 bp (phase 0) and [320,361) 41 bp (phase 2), 81 bp.
/// - txMirM (-) exons [100,160) [260,300) [350,380), thick [119,340) -> coding
///   pieces [260,300) 40 bp (phase 0) and [119,160) 41 bp (phase 2), 81 bp.
///   Transcript order is descending genomic, so [260,300) is its first coding
///   piece and [119,160) its last.
///
/// `gtf22` (GTF): the CDS rows lose the last three coding bases in transcript
/// order - [358,361) on +, [119,122) on - and 40 + 38 = 78 coding bases precede
/// the stop codon on both strands, so 78 % 3 = 0 and each codon is one row with
/// frame 0. UTR rows are the exon pieces outside the thick region, 5'/3' by
/// strand and ascending: [100,130) and [361,380) on +; [100,119) (3') and
/// [350,380) (5') on -.
///
/// GFF3 keeps 0.0.16's shape (`gtf22` is GTF-only): the CDS keeps the stop codon
/// and the codon rows carry no frame. With no map the gene id is the synthesized
/// `gene-<transcript id>` one.
#[test]
fn test_release_options_on_mirrored_pair() {
    let records = read_bed12("tests/data/bed12_mirror.bed", 1);
    assert_eq!(records.len(), 2, "mirrored pair: two records");
    let release = GxfOptions {
        source: Some("aprime"),
        utr_rows: true,
        gtf22: true,
        gff_escape: true,
        ..Default::default()
    };

    assert_eq!(
        as_text(records[0].to_gxf_with_options::<Gtf>(&release, None)),
        vec![
            "chrM1\taprime\tgene\t101\t380\t.\t+\t.\tgene_id \"txMirP\";",
            "chrM1\taprime\ttranscript\t101\t380\t.\t+\t.\tgene_id \"txMirP\"; transcript_id \"txMirP\";",
            "chrM1\taprime\texon\t101\t130\t.\t+\t.\tgene_id \"txMirP\"; transcript_id \"txMirP\"; exon_number \"1\";",
            "chrM1\taprime\texon\t181\t220\t.\t+\t.\tgene_id \"txMirP\"; transcript_id \"txMirP\"; exon_number \"2\";",
            "chrM1\taprime\texon\t321\t380\t.\t+\t.\tgene_id \"txMirP\"; transcript_id \"txMirP\"; exon_number \"3\";",
            "chrM1\taprime\tCDS\t181\t220\t.\t+\t0\tgene_id \"txMirP\"; transcript_id \"txMirP\"; exon_number \"2\";",
            "chrM1\taprime\tCDS\t321\t358\t.\t+\t2\tgene_id \"txMirP\"; transcript_id \"txMirP\"; exon_number \"3\";",
            "chrM1\taprime\tstart_codon\t181\t183\t.\t+\t0\tgene_id \"txMirP\"; transcript_id \"txMirP\"; exon_number \"2\";",
            "chrM1\taprime\tstop_codon\t359\t361\t.\t+\t0\tgene_id \"txMirP\"; transcript_id \"txMirP\"; exon_number \"3\";",
            "chrM1\taprime\tfive_prime_utr\t101\t130\t.\t+\t.\tgene_id \"txMirP\"; transcript_id \"txMirP\"; exon_number \"1\";",
            "chrM1\taprime\tthree_prime_utr\t362\t380\t.\t+\t.\tgene_id \"txMirP\"; transcript_id \"txMirP\"; exon_number \"3\";",
        ]
    );
    assert_eq!(
        as_text(records[1].to_gxf_with_options::<Gtf>(&release, None)),
        vec![
            "chrM1\taprime\tgene\t101\t380\t.\t-\t.\tgene_id \"txMirM\";",
            "chrM1\taprime\ttranscript\t101\t380\t.\t-\t.\tgene_id \"txMirM\"; transcript_id \"txMirM\";",
            "chrM1\taprime\texon\t101\t160\t.\t-\t.\tgene_id \"txMirM\"; transcript_id \"txMirM\"; exon_number \"3\";",
            "chrM1\taprime\texon\t261\t300\t.\t-\t.\tgene_id \"txMirM\"; transcript_id \"txMirM\"; exon_number \"2\";",
            "chrM1\taprime\texon\t351\t380\t.\t-\t.\tgene_id \"txMirM\"; transcript_id \"txMirM\"; exon_number \"1\";",
            "chrM1\taprime\tCDS\t123\t160\t.\t-\t2\tgene_id \"txMirM\"; transcript_id \"txMirM\"; exon_number \"3\";",
            "chrM1\taprime\tCDS\t261\t300\t.\t-\t0\tgene_id \"txMirM\"; transcript_id \"txMirM\"; exon_number \"2\";",
            "chrM1\taprime\tstart_codon\t298\t300\t.\t-\t0\tgene_id \"txMirM\"; transcript_id \"txMirM\"; exon_number \"2\";",
            "chrM1\taprime\tstop_codon\t120\t122\t.\t-\t0\tgene_id \"txMirM\"; transcript_id \"txMirM\"; exon_number \"3\";",
            "chrM1\taprime\tthree_prime_utr\t101\t119\t.\t-\t.\tgene_id \"txMirM\"; transcript_id \"txMirM\"; exon_number \"3\";",
            "chrM1\taprime\tfive_prime_utr\t351\t380\t.\t-\t.\tgene_id \"txMirM\"; transcript_id \"txMirM\"; exon_number \"1\";",
        ]
    );

    assert_eq!(
        as_text(records[0].to_gxf_with_options::<Gff>(&release, None)),
        vec![
            "chrM1\taprime\tgene\t101\t380\t.\t+\t.\tID=gene-txMirP;",
            "chrM1\taprime\tmRNA\t101\t380\t.\t+\t.\tID=txMirP;Parent=gene-txMirP;",
            "chrM1\taprime\texon\t101\t130\t.\t+\t.\tParent=txMirP;exon_number=1;",
            "chrM1\taprime\texon\t181\t220\t.\t+\t.\tParent=txMirP;exon_number=2;",
            "chrM1\taprime\texon\t321\t380\t.\t+\t.\tParent=txMirP;exon_number=3;",
            "chrM1\taprime\tCDS\t181\t220\t.\t+\t0\tParent=txMirP;exon_number=2;",
            "chrM1\taprime\tCDS\t321\t361\t.\t+\t2\tParent=txMirP;exon_number=3;",
            "chrM1\taprime\tstart_codon\t181\t183\t.\t+\t.\tParent=txMirP;exon_number=2;",
            "chrM1\taprime\tstop_codon\t359\t361\t.\t+\t.\tParent=txMirP;exon_number=3;",
            "chrM1\taprime\tfive_prime_UTR\t101\t130\t.\t+\t.\tParent=txMirP;exon_number=1;",
            "chrM1\taprime\tthree_prime_UTR\t362\t380\t.\t+\t.\tParent=txMirP;exon_number=3;",
        ]
    );
    assert_eq!(
        as_text(records[1].to_gxf_with_options::<Gff>(&release, None)),
        vec![
            "chrM1\taprime\tgene\t101\t380\t.\t-\t.\tID=gene-txMirM;",
            "chrM1\taprime\tmRNA\t101\t380\t.\t-\t.\tID=txMirM;Parent=gene-txMirM;",
            "chrM1\taprime\texon\t101\t160\t.\t-\t.\tParent=txMirM;exon_number=3;",
            "chrM1\taprime\texon\t261\t300\t.\t-\t.\tParent=txMirM;exon_number=2;",
            "chrM1\taprime\texon\t351\t380\t.\t-\t.\tParent=txMirM;exon_number=1;",
            "chrM1\taprime\tCDS\t120\t160\t.\t-\t2\tParent=txMirM;exon_number=3;",
            "chrM1\taprime\tCDS\t261\t300\t.\t-\t0\tParent=txMirM;exon_number=2;",
            "chrM1\taprime\tstart_codon\t298\t300\t.\t-\t.\tParent=txMirM;exon_number=2;",
            "chrM1\taprime\tstop_codon\t120\t122\t.\t-\t.\tParent=txMirM;exon_number=3;",
            "chrM1\taprime\tthree_prime_UTR\t101\t119\t.\t-\t.\tParent=txMirM;exon_number=3;",
            "chrM1\taprime\tfive_prime_UTR\t351\t380\t.\t-\t.\tParent=txMirM;exon_number=1;",
        ]
    );

    // GTF codon frame 0 on both strands: the codon is contiguous, so one row.
    for record in &records {
        for line in as_text(record.to_gxf_with_options::<Gtf>(&release, None)) {
            let columns: Vec<&str> = line.split('\t').collect();
            if columns[2] == "start_codon" || columns[2] == "stop_codon" {
                assert_eq!(columns[7], "0", "unframed codon row: {line}");
            }
        }
    }
}

/// A split stop codon on the - strand (`tests/data/bed12_split_stop.bed`), under
/// `gtf22`.
///
/// Derivation (0-based, half-open): exons [500,502) [570,598) [600,630), thick
/// [500,630) -> coding pieces [600,630) 30 bp, [570,598) 28 bp, [500,502) 2 bp
/// (60 bp, 20 codons). Transcript order is descending genomic, so the two last
/// coding bases are [500,502) and the third is the bottom base of [570,598):
/// 0.0.16 truncates the codon to that 2-base piece, while `gtf22` emits both
/// [570,571) and [500,502) as two `stop_codon` rows and drops [500,502) from the
/// CDS (it is emptied). Frames follow the phase rule: 57 coding bases precede
/// the [570,571) row (frame 0) and 58 precede [500,502) (frame 2), so the codon
/// spans frame 0 + frame 2.
#[test]
fn test_gtf22_split_stop_codon_on_minus() {
    let records = read_bed12("tests/data/bed12_split_stop.bed", 0);
    assert_eq!(records.len(), 1, "split stop codon: one record");
    let opts = GxfOptions {
        gtf22: true,
        ..Default::default()
    };

    assert_eq!(
        as_text(records[0].to_gxf_with_options::<Gtf>(&opts, None)),
        vec![
            "chrS\tgenepred\tgene\t501\t700\t.\t-\t.\tgene_id \"txSplit\";",
            "chrS\tgenepred\ttranscript\t501\t700\t.\t-\t.\tgene_id \"txSplit\"; transcript_id \"txSplit\";",
            "chrS\tgenepred\texon\t501\t502\t.\t-\t.\tgene_id \"txSplit\"; transcript_id \"txSplit\"; exon_number \"3\";",
            "chrS\tgenepred\texon\t571\t598\t.\t-\t.\tgene_id \"txSplit\"; transcript_id \"txSplit\"; exon_number \"2\";",
            "chrS\tgenepred\texon\t601\t630\t.\t-\t.\tgene_id \"txSplit\"; transcript_id \"txSplit\"; exon_number \"1\";",
            "chrS\tgenepred\tCDS\t572\t598\t.\t-\t0\tgene_id \"txSplit\"; transcript_id \"txSplit\"; exon_number \"2\";",
            "chrS\tgenepred\tCDS\t601\t630\t.\t-\t0\tgene_id \"txSplit\"; transcript_id \"txSplit\"; exon_number \"1\";",
            "chrS\tgenepred\tstart_codon\t628\t630\t.\t-\t0\tgene_id \"txSplit\"; transcript_id \"txSplit\"; exon_number \"1\";",
            "chrS\tgenepred\tstop_codon\t501\t502\t.\t-\t2\tgene_id \"txSplit\"; transcript_id \"txSplit\"; exon_number \"3\";",
            "chrS\tgenepred\tstop_codon\t571\t571\t.\t-\t0\tgene_id \"txSplit\"; transcript_id \"txSplit\"; exon_number \"2\";",
        ]
    );

    // Without `gtf22` (i.e. 0.0.16) the codon is truncated to the 2-base piece
    // and the CDS keeps all 60 coding bases.
    let text = as_text(records[0].to_gxf_with_options::<Gtf>(&GxfOptions::default(), None));
    assert!(text.contains(
        &"chrS\tgenepred\tstop_codon\t501\t502\t.\t-\t.\tgene_id \"txSplit\"; transcript_id \"txSplit\"; exon_number \"3\";".to_string()
    ));
    assert_eq!(
        text.iter()
            .filter(|line| line.split('\t').nth(2) == Some("CDS"))
            .count(),
        3,
        "0.0.16 keeps the emptied-by-gtf22 CDS piece"
    );
}

/// REVM1 MUST-FIX 1 and 2 on the review's own ± counterexamples
/// (`review017/data/review_cases.bed`, lines 1-2 and 5-6), under the release
/// options: a transcript with no coding base inside its exons, or with fewer
/// than three of them, is written as non-coding - gene / transcript / exon rows
/// only, no CDS, codon or UTR rows.
///
/// Derivation (0-based, half-open; the CDS is the exon ∩ [thickStart, thickEnd)
/// intersection):
///
/// * utrP / utrM (chrU 100-300, thick [180,220)): exons [100,120) and [250,270)
///   on +, [130,150) and [280,300) on -, so [180,220) is wholly intronic on
///   both strands and every coding piece is empty (0 coding bases). MUST-FIX 1
///   emitted one `five_prime_utr` and one `three_prime_utr` row per exon here.
/// * shortP / shortM (chrC 100-200, 2 bp CDS): one exon [130,132) on +, [168,170)
///   on -, i.e. 2 coding bases - fewer than the three a codon needs. MUST-FIX 2
///   emitted a 2 bp `start_codon` and a 2 bp `stop_codon` row here.
///
/// Three coding bases are a complete codon and are still written, and
/// `GxfOptions::default()` (no release option) keeps 0.0.16's rows for both
/// inputs, which is what the golden test pins.
#[test]
fn test_release_options_write_non_coding_without_exonic_cds_or_full_codon() {
    let release = GxfOptions {
        source: Some("aprime"),
        utr_rows: true,
        gtf22: true,
        gff_escape: true,
        ..Default::default()
    };

    // review_cases.bed line 1: chrU 100 300 utrP 0 + 180 220 0,0,0 2 20,20 0,150
    let utr_p = bed12_record(
        b"chrU",
        100,
        300,
        b"utrP",
        Strand::Forward,
        180,
        220,
        &[20, 20],
        &[0, 150],
    );
    // line 2, the same thick interval mirrored onto -: 2 20,20 30,180.
    let utr_m = bed12_record(
        b"chrU",
        100,
        300,
        b"utrM",
        Strand::Reverse,
        180,
        220,
        &[20, 20],
        &[30, 180],
    );
    // line 5: chrC 100 200 shortP 0 + 130 132 0,0,0 1 2 30
    let short_p = bed12_record(
        b"chrC",
        100,
        200,
        b"shortP",
        Strand::Forward,
        130,
        132,
        &[2],
        &[30],
    );
    // line 6, the same mirrored onto -: thick 168,170, exon 2 at 68.
    let short_m = bed12_record(
        b"chrC",
        100,
        200,
        b"shortM",
        Strand::Reverse,
        168,
        170,
        &[2],
        &[68],
    );

    assert_eq!(
        as_text(utr_p.to_gxf_with_options::<Gtf>(&release, None)),
        vec![
            "chrU\taprime\tgene\t101\t300\t.\t+\t.\tgene_id \"utrP\";",
            "chrU\taprime\ttranscript\t101\t300\t.\t+\t.\tgene_id \"utrP\"; transcript_id \"utrP\";",
            "chrU\taprime\texon\t101\t120\t.\t+\t.\tgene_id \"utrP\"; transcript_id \"utrP\"; exon_number \"1\";",
            "chrU\taprime\texon\t251\t270\t.\t+\t.\tgene_id \"utrP\"; transcript_id \"utrP\"; exon_number \"2\";",
        ]
    );
    assert_eq!(
        as_text(utr_m.to_gxf_with_options::<Gtf>(&release, None)),
        vec![
            "chrU\taprime\tgene\t101\t300\t.\t-\t.\tgene_id \"utrM\";",
            "chrU\taprime\ttranscript\t101\t300\t.\t-\t.\tgene_id \"utrM\"; transcript_id \"utrM\";",
            "chrU\taprime\texon\t131\t150\t.\t-\t.\tgene_id \"utrM\"; transcript_id \"utrM\"; exon_number \"2\";",
            "chrU\taprime\texon\t281\t300\t.\t-\t.\tgene_id \"utrM\"; transcript_id \"utrM\"; exon_number \"1\";",
        ]
    );
    assert_eq!(
        as_text(short_p.to_gxf_with_options::<Gtf>(&release, None)),
        vec![
            "chrC\taprime\tgene\t101\t200\t.\t+\t.\tgene_id \"shortP\";",
            "chrC\taprime\ttranscript\t101\t200\t.\t+\t.\tgene_id \"shortP\"; transcript_id \"shortP\";",
            "chrC\taprime\texon\t131\t132\t.\t+\t.\tgene_id \"shortP\"; transcript_id \"shortP\"; exon_number \"1\";",
        ]
    );
    assert_eq!(
        as_text(short_m.to_gxf_with_options::<Gtf>(&release, None)),
        vec![
            "chrC\taprime\tgene\t101\t200\t.\t-\t.\tgene_id \"shortM\";",
            "chrC\taprime\ttranscript\t101\t200\t.\t-\t.\tgene_id \"shortM\"; transcript_id \"shortM\";",
            "chrC\taprime\texon\t169\t170\t.\t-\t.\tgene_id \"shortM\"; transcript_id \"shortM\"; exon_number \"1\";",
        ]
    );

    // The gate belongs to each release option on its own, so a caller setting
    // only one of them never gets a UTR, CDS or codon row for these records.
    for options in [
        release,
        GxfOptions {
            utr_rows: true,
            gtf22: false,
            ..release
        },
        GxfOptions {
            utr_rows: false,
            gtf22: true,
            ..release
        },
    ] {
        for record in [&utr_p, &utr_m, &short_p, &short_m] {
            for line in as_text(record.to_gxf_with_options::<Gtf>(&options, None)) {
                let feature = line.split('\t').nth(2).unwrap();
                assert!(
                    !matches!(
                        feature,
                        "CDS" | "start_codon" | "stop_codon" | "five_prime_utr" | "three_prime_utr"
                    ),
                    "non-coding transcript wrote a {feature} row: {line}"
                );
            }
            // GFF3 is gated the same way, only the UTR feature names differ.
            for line in as_text(record.to_gxf_with_options::<Gff>(&options, None)) {
                let feature = line.split('\t').nth(2).unwrap();
                assert!(
                    !matches!(
                        feature,
                        "CDS" | "start_codon" | "stop_codon" | "five_prime_UTR" | "three_prime_UTR"
                    ),
                    "non-coding transcript wrote a {feature} row: {line}"
                );
            }
        }
    }

    // Three coding bases are one complete codon, so they are still written (and
    // `gtf22` trims all three out of the CDS rows into the codon rows).
    let three_p = bed12_record(
        b"chrC",
        100,
        200,
        b"threeP",
        Strand::Forward,
        130,
        133,
        &[3],
        &[30],
    );
    let three_rows = as_text(three_p.to_gxf_with_options::<Gtf>(&release, None));
    assert!(three_rows.contains(
        &"chrC\taprime\tstart_codon\t131\t133\t.\t+\t0\tgene_id \"threeP\"; transcript_id \"threeP\"; exon_number \"1\";".to_string()
    ));
    assert!(three_rows.contains(
        &"chrC\taprime\tstop_codon\t131\t133\t.\t+\t0\tgene_id \"threeP\"; transcript_id \"threeP\"; exon_number \"1\";".to_string()
    ));
    assert_eq!(
        three_rows
            .iter()
            .filter(|line| line.split('\t').nth(2) == Some("CDS"))
            .count(),
        0,
        "the 3 bp CDS is trimmed into the codon rows, not dropped as truncated"
    );

    // Default options are not gated: 0.0.16 keeps the 2 bp CDS row and its
    // truncated codon rows, which is the golden test's behaviour.
    let legacy = as_text(short_p.to_gxf_with_options::<Gtf>(&GxfOptions::default(), None));
    assert!(legacy.contains(
        &"chrC\tgenepred\tCDS\t131\t132\t.\t+\t0\tgene_id \"shortP\"; transcript_id \"shortP\"; exon_number \"1\";".to_string()
    ));
    assert!(legacy.contains(
        &"chrC\tgenepred\tstart_codon\t131\t132\t.\t+\t.\tgene_id \"shortP\"; transcript_id \"shortP\"; exon_number \"1\";".to_string()
    ));
    assert!(legacy.contains(
        &"chrC\tgenepred\tstop_codon\t131\t132\t.\t+\t.\tgene_id \"shortP\"; transcript_id \"shortP\"; exon_number \"1\";".to_string()
    ));
}

/// Builds a `GenePred` from the BED12 fields of one review counterexample line
/// (0-based, half-open), the way `read_bed12` reads them.
#[allow(clippy::too_many_arguments)]
fn bed12_record(
    chrom: &[u8],
    start: u64,
    end: u64,
    name: &[u8],
    strand: Strand,
    thick_start: u64,
    thick_end: u64,
    sizes: &[u32],
    starts: &[u32],
) -> GenePred {
    Bed12 {
        chrom: chrom.to_vec(),
        start,
        end,
        name: name.to_vec(),
        score: 0,
        strand,
        thick_start,
        thick_end,
        item_rgb: genepred::bed::Rgb(0, 0, 0),
        block_count: sizes.len() as u32,
        block_sizes: sizes.to_vec(),
        block_starts: starts.to_vec(),
        extras: Extras::new(),
    }
    .into()
}

/// SO percent-encoding of a seqid and its IDs (`tests/data/bed12_special_seqid.bed`,
/// whose seqid and transcript name are `c;1=a,b&%x#y`), GFF3 only.
///
/// The two rules differ by design (release_design.md sec.3): the seqid escapes
/// every byte outside `[a-zA-Z0-9.:^*$@!+_?-|]`, while a column-9 value only
/// escapes `%`, `;`, `=`, `&`, `,` and the control bytes. The only byte of this
/// name the rules disagree on is `#`, which is fine inside column 9 - only a
/// `#` starting a line (or column 1) opens a GFF3 comment - so it stays literal
/// in the ID and GTF writes the name verbatim throughout.
#[test]
fn test_gff_escape_seqid_and_id() {
    let records = read_bed12("tests/data/bed12_special_seqid.bed", 0);
    assert_eq!(records.len(), 1, "special seqid: one record");
    let opts = GxfOptions {
        gff_escape: true,
        ..Default::default()
    };

    let gff = as_text(records[0].to_gxf_with_options::<Gff>(&opts, None));
    let seqid = "c%3B1%3Da%2Cb%26%25x%23y";
    assert_eq!(
        gff[0],
        format!("{seqid}\tgenepred\tgene\t11\t210\t.\t+\t.\tID=gene-c%3B1%3Da%2Cb%26%25x#y;")
    );
    assert_eq!(
        gff[1],
        format!("{seqid}\tgenepred\tmRNA\t11\t210\t.\t+\t.\tID=c%3B1%3Da%2Cb%26%25x#y;Parent=gene-c%3B1%3Da%2Cb%26%25x#y;")
    );
    assert!(
        gff.iter()
            .all(|line| line.starts_with(&format!("{seqid}\t"))),
        "every row carries the escaped seqid"
    );

    let gtf = as_text(records[0].to_gxf_with_options::<Gtf>(&opts, None));
    assert!(
        gtf.iter().all(|line| line.starts_with("c;1=a,b&%x#y\t")),
        "GTF has no escaping, even with gff_escape set"
    );
}

/// STEP 1 of M1 ("Default byte identity"): every GXF path the release design
/// lists, replayed over every golden input with `GxfOptions::default()`, must
/// rebuild `tests/data/golden_0_0_16.{gtf,gff}` — the frozen output of the
/// published `genepred = "=0.0.16"` (registry checksum recorded in the golden
/// commit's message) — byte for byte. Any of the four new options that moves a
/// default byte fails here.
#[test]
fn test_gxf_options_default_rebuilds_golden_0_0_16_bytes() {
    let (gtf, gff) = replay_golden_sequence();

    let expected_gtf: &[u8] = include_bytes!("data/golden_0_0_16.gtf");
    let expected_gff: &[u8] = include_bytes!("data/golden_0_0_16.gff");
    assert_eq!(gtf.as_slice(), expected_gtf, "GTF golden mismatch");
    assert_eq!(gff.as_slice(), expected_gff, "GFF golden mismatch");
}

/// Golden inputs, in the contract's fixed order, with the number of trailing BED
/// extras read into each record (extra keys `13`, `14`, ...).
const GOLDEN_INPUTS: &[(&str, usize)] = &[
    ("bed12.bed", 0),
    ("bed12_extra.bed", 2),
    ("bed12_mirror.bed", 1),
    ("bed12_split_stop.bed", 0),
    ("bed12_special_seqid.bed", 0),
];

/// Number of GXF paths replayed per record per map variant: `to_gxf`,
/// `to_gxf_with_additional_fields` with 0 and with 1, the four
/// `to_gxf_with_options` combinations, and `to_gxf_gene_line`.
const GOLDEN_PATHS: usize = 8;

/// Rebuilds the frozen golden sequence: for each input, for each GXF path, with
/// and then without the transcript -> gene map, plus the `Writer` path over the
/// whole input.
fn replay_golden_sequence() -> (Vec<u8>, Vec<u8>) {
    let mut gtf = Vec::new();
    let mut gff = Vec::new();

    for (name, additional_fields) in GOLDEN_INPUTS {
        let records = read_bed12(&format!("tests/data/{name}"), *additional_fields);
        assert!(!records.is_empty(), "golden input has no records: {name}");
        let map = golden_gene_map(&records);
        let matrix = [
            GeneLine::Include,
            GeneLine::Include,
            GeneLine::Omit,
            GeneLine::Omit,
        ]
        .into_iter()
        .zip([
            TranscriptParent::Include,
            TranscriptParent::Omit,
            TranscriptParent::Include,
            TranscriptParent::Omit,
        ])
        .map(|(gene_line, transcript_parent)| GxfOptions {
            gene_line,
            transcript_parent,
            ..Default::default()
        })
        .collect::<Vec<_>>();

        for with_map in [true, false] {
            let map = with_map.then_some(&map);
            for record in &records {
                for path in 0..GOLDEN_PATHS {
                    replay_golden_path::<Gtf>(record, path, &matrix, map, &mut gtf);
                    replay_golden_path::<Gff>(record, path, &matrix, map, &mut gff);
                }
            }

            let mut buffer = Vec::new();
            Writer::<Gtf>::from_records_with_options(
                &records,
                &mut buffer,
                &WriterOptions::default(),
            )
            .unwrap();
            gtf.extend_from_slice(&buffer);
            let mut buffer = Vec::new();
            Writer::<Gff>::from_records_with_options(
                &records,
                &mut buffer,
                &WriterOptions::default(),
            )
            .unwrap();
            gff.extend_from_slice(&buffer);
        }
    }

    (gtf, gff)
}

/// The transcript -> gene map the golden files were written with: each record's
/// name maps to `gene_<name>`, derived only from the input so both the writer
/// and this replay build the same map.
fn golden_gene_map(records: &[GenePred]) -> HashMap<String, String> {
    records
        .iter()
        .filter_map(|record| {
            let name = record
                .name()
                .map(|name| String::from_utf8_lossy(name).into_owned())?;
            Some((name.clone(), format!("gene_{name}")))
        })
        .collect()
}

/// Whether a record carries at least one numeric extra, i.e. whether
/// `to_gxf_with_additional_fields(1, _)` applies to it.
fn has_numeric_extra(record: &GenePred) -> bool {
    record
        .extras()
        .keys()
        .any(|key| std::str::from_utf8(key).is_ok_and(|key| key.parse::<u64>().is_ok()))
}

/// Appends one GXF path's rows to `out`.
fn replay_golden_path<K>(
    record: &GenePred,
    path: usize,
    matrix: &[GxfOptions],
    map: Option<&HashMap<String, String>>,
    out: &mut Vec<u8>,
) where
    K: BedFormat,
{
    let lines = match path {
        0 => record.to_gxf::<K>(map),
        1 => record.to_gxf_with_additional_fields::<K>(0, map),
        2 => {
            if has_numeric_extra(record) {
                record.to_gxf_with_additional_fields::<K>(1, map)
            } else {
                return;
            }
        }
        3..=6 => record.to_gxf_with_options::<K>(&matrix[path - 3], map),
        7 => vec![record.to_gxf_gene_line::<K>(0, map)],
        _ => unreachable!("golden GXF path"),
    };
    for line in lines {
        out.extend_from_slice(&line);
        out.push(b'\n');
    }
}

/// Reads a BED12 file with `additional_fields` trailing extras.
fn read_bed12(path: &str, additional_fields: usize) -> Vec<GenePred> {
    let options = ReaderOptions::new().additional_fields(additional_fields);
    let mut reader = Reader::<Bed12>::from_path_with_custom_fields(path, options)
        .unwrap_or_else(|error| panic!("reading {path}: {error}"));
    reader
        .records()
        .map(|record| record.unwrap_or_else(|error| panic!("record in {path}: {error}")))
        .collect()
}
