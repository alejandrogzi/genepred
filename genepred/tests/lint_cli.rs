use std::{fs, process::Command};

use tempfile::tempdir;

/// Returns the path to the compiled `genepred` test binary.
fn genepred() -> &'static str {
    env!("CARGO_BIN_EXE_genepred")
}

/// Verifies valid BED input exits successfully without diagnostics.
#[test]
fn lint_valid_bed_exits_zero() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("valid.bed");
    fs::write(&path, "chr1\t0\t10\nchr1\t20\t30\n").unwrap();

    let output = Command::new(genepred())
        .args(["lint", path.to_str().unwrap()])
        .output()
        .unwrap();

    assert!(output.status.success());
}

/// Verifies invalid BED input exits with code 1.
#[test]
fn lint_invalid_bed_exits_one() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("invalid.bed");
    fs::write(&path, "chr1\t0\t10\nchr1\t20\t20\n").unwrap();

    let output = Command::new(genepred())
        .args(["lint", path.to_str().unwrap()])
        .output()
        .unwrap();

    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("invalid=1"));
}

/// Verifies warning mode reports diagnostics but exits successfully.
#[test]
fn lint_warn_invalid_bed_exits_zero() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("invalid.bed");
    fs::write(&path, "chr1\t0\t10\nchr1\t20\t20\n").unwrap();

    let output = Command::new(genepred())
        .args(["lint", "--warn", path.to_str().unwrap()])
        .output()
        .unwrap();

    assert!(output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("WARN:"));
    assert!(stderr.contains("invalid=1"));
}

/// Verifies BED pruning emits only valid original records.
#[test]
fn lint_prune_bed_outputs_only_valid_original_records() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("mixed.bed");
    fs::write(&path, "# comment\nchr1\t0\t10\nchr1\t20\t20\nchr2\t5\t8\n").unwrap();

    let output = Command::new(genepred())
        .args(["lint", "--prune", path.to_str().unwrap()])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "chr1\t0\t10\nchr2\t5\t8\n"
    );
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("PRUNE:"));
}

/// Verifies `--additional-fields` resolves an ambiguous BED6 layout.
#[test]
fn lint_additional_fields_disambiguates_bed6_extra_columns() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("bed6_plus4.bed");
    fs::write(&path, "chr1\t0\t10\ttx1\t500\t+\t0\t10\tnot_rgb\textra\n").unwrap();

    let output = Command::new(genepred())
        .args(["lint", "-a", "4", path.to_str().unwrap()])
        .output()
        .unwrap();

    assert!(output.status.success());
}

/// Verifies BED pruning preserves original lines with additional fields.
#[test]
fn lint_prune_additional_fields_preserves_original_bed_line() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("bed6_plus4.bed");
    let content = "chr1\t0\t10\ttx1\t500\t+\t0\t10\tnot_rgb\textra\n";
    fs::write(&path, content).unwrap();

    let output = Command::new(genepred())
        .args([
            "lint",
            "--prune",
            "--additional-fields",
            "4",
            path.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert_eq!(String::from_utf8(output.stdout).unwrap(), content);
}

/// Verifies BED additional fields are rejected for GTF input.
#[test]
fn lint_rejects_additional_fields_for_gtf() {
    let output = Command::new(genepred())
        .args(["lint", "-a", "1", "tests/data/simple.gtf"])
        .output()
        .unwrap();

    assert_eq!(output.status.code(), Some(2));
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("--additional-fields/-a is only supported for BED input"));
    assert!(stderr.contains("detected GTF"));
}

/// Verifies BED additional fields are rejected for GFF prune mode.
#[test]
fn lint_prune_rejects_additional_fields_for_gff() {
    let output = Command::new(genepred())
        .args(["lint", "--prune", "-a", "1", "tests/data/simple.gff"])
        .output()
        .unwrap();

    assert_eq!(output.status.code(), Some(2));
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("--additional-fields/-a is only supported for BED input"));
    assert!(stderr.contains("detected GFF"));
}

/// Verifies invalid BED9 RGB values are lint diagnostics.
#[test]
fn lint_invalid_bed9_rgb_exits_one() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("invalid_rgb.bed");
    fs::write(&path, "chr1\t0\t10\ttx1\t500\t+\t0\t10\t255,0,300\n").unwrap();

    let output = Command::new(genepred())
        .args(["lint", path.to_str().unwrap()])
        .output()
        .unwrap();

    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("itemRgb"));
    assert!(stderr.contains("invalid=1"));
}

/// Verifies invalid BED12 RGB values are lint diagnostics.
#[test]
fn lint_invalid_bed12_rgb_exits_one() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("invalid_rgb.bed");
    fs::write(
        &path,
        "chr1\t0\t10\ttx1\t500\t+\t0\t10\t255,0,bad\t1\t10\t0\n",
    )
    .unwrap();

    let output = Command::new(genepred())
        .args(["lint", path.to_str().unwrap()])
        .output()
        .unwrap();

    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("itemRgb"));
    assert!(stderr.contains("invalid=1"));
}

/// Verifies unsupported explicit BED layouts are rejected.
#[test]
fn lint_rejects_unsupported_explicit_bed_layout() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("bed10.bed");
    fs::write(&path, "chr1\t0\t10\ttx1\t500\t+\t0\t10\t0,0,0\textra\n").unwrap();

    let output = Command::new(genepred())
        .args(["lint", "-a", "3", path.to_str().unwrap()])
        .output()
        .unwrap();

    assert_eq!(output.status.code(), Some(2));
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("leaves 7 standard field(s)"));
}

/// Verifies `--warn` and `--prune` are mutually exclusive.
#[test]
fn lint_warn_and_prune_conflict() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("valid.bed");
    fs::write(&path, "chr1\t0\t10\n").unwrap();

    let output = Command::new(genepred())
        .args(["lint", "--warn", "--prune", path.to_str().unwrap()])
        .output()
        .unwrap();

    assert_eq!(output.status.code(), Some(2));
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("cannot be used with"));
}

/// Verifies gzip BED input uses the buffered fallback path.
#[cfg(feature = "gzip")]
#[test]
fn lint_gzip_bed_uses_buffered_fallback() {
    let output = Command::new(genepred())
        .args(["lint", "tests/data/bed3.bed.gz"])
        .output()
        .unwrap();

    assert!(output.status.success());
}

/// Verifies gzip GTF prune uses the spooled fallback path.
#[cfg(feature = "gzip")]
#[test]
fn lint_prune_gzip_gtf_uses_spooled_fallback() {
    let output = Command::new(genepred())
        .args(["lint", "--prune", "tests/data/simple.gtf.gz"])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        fs::read_to_string("tests/data/simple.gtf").unwrap()
    );
}

/// Verifies GFF prune preserves comments and valid transcript features.
#[test]
fn lint_prune_gff_preserves_comments_and_valid_transcripts() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("mixed.gff3");
    fs::write(
        &path,
        concat!(
            "##gff-version 3\n",
            "# keep this\n",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Name=Tx1\n",
            "chr1\tsrc\texon\t1\t50\t.\t+\t.\tParent=tx1\n",
            "chr1\tsrc\tmRNA\t200\t250\t.\t+\t.\tID=bad;Name=Bad\n",
            "chr1\tsrc\texon\t100\t120\t.\t+\t.\tParent=bad\n",
        ),
    )
    .unwrap();

    let output = Command::new(genepred())
        .args(["lint", "--prune", path.to_str().unwrap()])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        concat!(
            "##gff-version 3\n",
            "# keep this\n",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Name=Tx1\n",
            "chr1\tsrc\texon\t1\t50\t.\t+\t.\tParent=tx1\n",
        )
    );
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("PRUNE:"));
    assert!(stderr.contains("invalid=1"));
}

/// Verifies GTF prune preserves comments and valid transcript features.
#[test]
fn lint_prune_gtf_preserves_comments_and_valid_transcripts() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("mixed.gtf");
    fs::write(
        &path,
        concat!(
            "# keep this\n",
            "chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            "chr1\tsrc\texon\t1\t50\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            "chr1\tsrc\ttranscript\t200\t250\t.\t+\t.\tgene_id \"g2\"; transcript_id \"bad\";\n",
            "chr1\tsrc\texon\t100\t120\t.\t+\t.\tgene_id \"g2\"; transcript_id \"bad\";\n",
        ),
    )
    .unwrap();

    let output = Command::new(genepred())
        .args(["lint", "--prune", path.to_str().unwrap()])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        concat!(
            "# keep this\n",
            "chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            "chr1\tsrc\texon\t1\t50\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
        )
    );
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("PRUNE:"));
    assert!(stderr.contains("invalid=1"));
}
