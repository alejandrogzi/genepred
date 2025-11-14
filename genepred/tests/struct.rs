use genepred::bed::{Bed12, Bed3, Bed4, Bed6, Rgb, Strand};
use genepred::GenePred;

#[test]
fn test_genepred_from_coords() {
    let gene = GenePred::from_coords(
        "chr1".to_string(),
        10,
        20,
        vec!["extra1".to_string(), "extra2".to_string()],
    );
    assert_eq!(gene.chrom(), "chr1");
    assert_eq!(gene.start(), 10);
    assert_eq!(gene.end(), 20);
    assert_eq!(
        gene.extras(),
        &vec!["extra1".to_string(), "extra2".to_string()]
    );
    assert!(gene.name().is_none());
    assert!(gene.score().is_none());
}

#[test]
fn test_genepred_from_bed3() {
    let bed3 = Bed3 {
        chrom: "chr1".to_string(),
        start: 10,
        end: 20,
        extras: vec![],
    };
    let gene: GenePred = bed3.into();
    assert_eq!(gene.chrom(), "chr1");
    assert_eq!(gene.start(), 10);
    assert_eq!(gene.end(), 20);
    assert!(gene.name().is_none());
}

#[test]
fn test_genepred_from_bed4() {
    let bed4 = Bed4 {
        chrom: "chr1".to_string(),
        start: 10,
        end: 20,
        name: "geneA".to_string(),
        extras: vec![],
    };
    let gene: GenePred = bed4.into();
    assert_eq!(gene.chrom(), "chr1");
    assert_eq!(gene.start(), 10);
    assert_eq!(gene.end(), 20);
    assert_eq!(gene.name().unwrap(), "geneA");
}

#[test]
fn test_genepred_from_bed6() {
    let bed6 = Bed6 {
        chrom: "chr1".to_string(),
        start: 10,
        end: 20,
        name: "geneA".to_string(),
        score: 100,
        strand: Strand::Forward,
        extras: vec![],
    };
    let gene: GenePred = bed6.into();
    assert_eq!(gene.chrom(), "chr1");
    assert_eq!(gene.start(), 10);
    assert_eq!(gene.end(), 20);
    assert_eq!(gene.name().unwrap(), "geneA");
    assert_eq!(gene.strand().unwrap(), Strand::Forward);
}

#[test]
fn test_genepred_from_bed12() {
    let bed12 = Bed12 {
        chrom: "chr1".to_string(),
        start: 10,
        end: 100,
        name: "geneA".to_string(),
        score: 1000,
        strand: Strand::Forward,
        thick_start: 10,
        thick_end: 100,
        item_rgb: Rgb(255, 0, 0),
        block_count: 2,
        block_sizes: vec![10, 20],
        block_starts: vec![0, 30],
        extras: vec![],
    };
    let gene: GenePred = bed12.into();
    assert_eq!(gene.chrom(), "chr1");
    assert_eq!(gene.start(), 10);
    assert_eq!(gene.end(), 100);
    assert_eq!(gene.name().unwrap(), "geneA");
    assert_eq!(gene.score().unwrap(), 1000);
    assert_eq!(gene.strand().unwrap(), Strand::Forward);
    assert_eq!(gene.thick_start().unwrap(), 10);
    assert_eq!(gene.thick_end().unwrap(), 100);
    assert_eq!(gene.item_rgb().unwrap(), Rgb(255, 0, 0));
    assert_eq!(gene.block_count().unwrap(), 2);
    assert_eq!(gene.block_sizes().unwrap(), &[10, 20]);
    assert_eq!(gene.block_starts().unwrap(), &[0, 30]);
}

#[test]
fn test_genepred_len_is_empty() {
    let gene = GenePred::from_coords("chr1".to_string(), 10, 20, vec![]);
    assert_eq!(gene.len(), 10);
    assert!(!gene.is_empty());

    let empty_gene = GenePred::from_coords("chr1".to_string(), 10, 10, vec![]);
    assert_eq!(empty_gene.len(), 0);
    assert!(empty_gene.is_empty());

    let inverted_gene = GenePred::from_coords("chr1".to_string(), 20, 10, vec![]);
    assert_eq!(inverted_gene.len(), 0);
    assert!(inverted_gene.is_empty());
}

#[test]
fn test_genepred_setters() {
    let mut gene = GenePred::from_coords("chr1".to_string(), 10, 20, vec![]);
    gene.set_chrom("chrX".to_string());
    gene.set_start(5);
    gene.set_end(25);
    gene.set_name(Some("new_gene".to_string()));
    gene.set_score(Some(500));
    gene.set_strand(Some(Strand::Reverse));
    gene.set_thick_start(Some(6));
    gene.set_thick_end(Some(24));
    gene.set_item_rgb(Some(Rgb(0, 0, 255)));
    gene.set_block_count(Some(1));
    gene.set_block_sizes(Some(vec![19]));
    gene.set_block_starts(Some(vec![1]));
    gene.set_extras(vec!["new_extra".to_string()]);

    assert_eq!(gene.chrom(), "chrX");
    assert_eq!(gene.start(), 5);
    assert_eq!(gene.end(), 25);
    assert_eq!(gene.name().unwrap(), "new_gene");
    assert_eq!(gene.score().unwrap(), 500);
    assert_eq!(gene.strand().unwrap(), Strand::Reverse);
    assert_eq!(gene.thick_start().unwrap(), 6);
    assert_eq!(gene.thick_end().unwrap(), 24);
    assert_eq!(gene.item_rgb().unwrap(), Rgb(0, 0, 255));
    assert_eq!(gene.block_count().unwrap(), 1);
    assert_eq!(gene.block_sizes().unwrap(), &[19]);
    assert_eq!(gene.block_starts().unwrap(), &[1]);
    assert_eq!(gene.extras(), &vec!["new_extra".to_string()]);

    gene.add_extra("another_extra".to_string());
    assert_eq!(
        gene.extras(),
        &vec!["new_extra".to_string(), "another_extra".to_string()]
    );
    gene.clear_extras();
    assert!(gene.extras().is_empty());
}

#[test]
fn test_genepred_exons() {
    // No blocks defined
    let gene1 = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    assert_eq!(gene1.exons(), vec![(10, 100)]);

    // With blocks defined
    let mut gene2 = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    gene2.set_block_count(Some(2));
    gene2.set_block_sizes(Some(vec![10, 20]));
    gene2.set_block_starts(Some(vec![0, 30])); // Relative to start (10)
    assert_eq!(gene2.exons(), vec![(10, 20), (40, 60)]);

    // With empty blocks
    let mut gene3 = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    gene3.set_block_count(Some(0));
    gene3.set_block_sizes(Some(vec![]));
    gene3.set_block_starts(Some(vec![]));
    assert_eq!(gene3.exons(), vec![(10, 100)]);
}

#[test]
fn test_genepred_introns() {
    // No introns (single exon)
    let gene1 = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    assert!(gene1.introns().is_empty());

    // With multiple exons
    let mut gene2 = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    gene2.set_block_count(Some(2));
    gene2.set_block_sizes(Some(vec![10, 20]));
    gene2.set_block_starts(Some(vec![0, 30])); // Exons: (10,20), (40,60)
    assert_eq!(gene2.introns(), vec![(20, 40)]);

    // With more exons
    let mut gene3 = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    gene3.set_block_count(Some(3));
    gene3.set_block_sizes(Some(vec![10, 10, 10]));
    gene3.set_block_starts(Some(vec![0, 20, 40])); // Exons: (10,20), (30,40), (50,60)
    assert_eq!(gene3.introns(), vec![(20, 30), (40, 50)]);
}

#[test]
fn test_genepred_exonic_intronic_length() {
    let mut gene = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    gene.set_block_count(Some(2));
    gene.set_block_sizes(Some(vec![10, 20]));
    gene.set_block_starts(Some(vec![0, 30])); // Exons: (10,20), (40,60)
                                              // Exonic length: (20-10) + (60-40) = 10 + 20 = 30
                                              // Intronic length: (40-20) = 20
    assert_eq!(gene.exonic_length(), 30);
    assert_eq!(gene.intronic_length(), 20);

    let gene_no_blocks = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    assert_eq!(gene_no_blocks.exonic_length(), 90); // (100-10)
    assert_eq!(gene_no_blocks.intronic_length(), 0);
}

#[test]
fn test_genepred_coding_exons_cds_length() {
    let mut gene = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    gene.set_block_count(Some(2));
    gene.set_block_sizes(Some(vec![10, 20]));
    gene.set_block_starts(Some(vec![0, 30])); // Exons: (10,20), (40,60)

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
    let mut gene = GenePred::from_coords("chr1".to_string(), 10, 20, vec![]);
    gene.add_extra("tag1=value1;tag2=value2".to_string());
    gene.add_extra("tag3=value3".to_string());

    assert_eq!(
        gene.unnest_extras(";"),
        vec![
            "tag1=value1".to_string(),
            "tag2=value2".to_string(),
            "tag3=value3".to_string()
        ]
    );

    assert_eq!(
        gene.unnest_extras("="),
        vec![
            "tag1".to_string(),
            "value1;tag2".to_string(),
            "value2".to_string(),
            "tag3".to_string(),
            "value3".to_string()
        ]
    );
}

#[test]
fn test_genepred_overlaps() {
    let gene = GenePred::from_coords("chr1".to_string(), 10, 20, vec![]);

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
    let mut gene = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    gene.set_block_count(Some(2));
    gene.set_block_sizes(Some(vec![10, 20]));
    gene.set_block_starts(Some(vec![0, 30])); // Exons: (10,20), (40,60)

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
    let gene1 = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    assert_eq!(gene1.exon_count(), 1);
    assert_eq!(gene1.intron_count(), 0);

    let mut gene2 = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    gene2.set_block_count(Some(2));
    gene2.set_block_sizes(Some(vec![10, 20]));
    gene2.set_block_starts(Some(vec![0, 30])); // Exons: (10,20), (40,60)
    assert_eq!(gene2.exon_count(), 2);
    assert_eq!(gene2.intron_count(), 1);

    let mut gene3 = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    gene3.set_block_count(Some(0));
    gene3.set_block_sizes(Some(vec![]));
    gene3.set_block_starts(Some(vec![]));
    assert_eq!(gene3.exon_count(), 1);
    assert_eq!(gene3.intron_count(), 0);
}

#[test]
fn test_genepred_display() {
    let gene1 = GenePred::from_coords("chr1".to_string(), 10, 20, vec![]);
    assert_eq!(format!("{}", gene1), "chr1\t10\t20");

    let mut gene2 = GenePred::from_coords("chr1".to_string(), 10, 20, vec![]);
    gene2.set_name(Some("geneA".to_string()));
    gene2.set_score(Some(100));
    gene2.set_strand(Some(Strand::Forward));
    assert_eq!(format!("{}", gene2), "chr1\t10\t20\tgeneA\t100\t+");

    let mut gene3 = GenePred::from_coords("chr1".to_string(), 10, 100, vec![]);
    gene3.set_name(Some("geneB".to_string()));
    gene3.set_score(Some(500));
    gene3.set_strand(Some(Strand::Reverse));
    gene3.set_thick_start(Some(10));
    gene3.set_thick_end(Some(100));
    gene3.set_item_rgb(Some(Rgb(255, 0, 0)));
    gene3.set_block_count(Some(2));
    gene3.set_block_sizes(Some(vec![10, 20]));
    gene3.set_block_starts(Some(vec![0, 30]));
    gene3.set_extras(vec!["extra1".to_string(), "extra2".to_string()]);
    assert_eq!(
        format!("{}", gene3),
        "chr1\t10\t100\tgeneB\t500\t-\t10\t100\t255,0,0\t2\t10,20\t0,30\textra1\textra2"
    );
}
