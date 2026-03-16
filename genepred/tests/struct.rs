use std::collections::HashMap;

use genepred::bed::{Bed12, Bed3, Bed4, Bed5, Bed6, Bed8, Bed9};
use genepred::{ExtraValue, Extras, GenePred, Gff, Gtf, Strand};

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
