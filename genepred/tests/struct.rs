use genepred::bed::{Bed12, Bed3, Bed4, Bed6};
use genepred::{ExtraValue, Extras, GenePred, Strand};

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
