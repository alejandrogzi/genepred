## Examples

### Reading a BED3 File

BED3 is the minimal format with chromosome, start, and end coordinates:

```rust,no_run
use genepred::{Reader, Bed3};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = Reader::<Bed3>::from_path("data/regions.bed")?;
    for record in reader.records() {
        let record = record?;
        println!("Region: {}:{}-{}",
            String::from_utf8_lossy(&record.chrom),
            record.start,
            record.end);
    }
    Ok(())
}
```

### Reading a BED6 File

BED6 includes name, score, and strand information:

```rust,no_run
use genepred::{Reader, Bed6};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = Reader::<Bed6>::from_path("data/genes.bed")?;
    for record in reader.records() {
        let record = record?;
        let name = record.name.unwrap_or(b"<unknown>".to_vec());
        let strand = record
            .strand
            .map(|s| s.to_string())
            .unwrap_or_else(|| ".".to_string());
        println!(
            "Gene: {} on {} strand at {}:{}-{}",
            String::from_utf8_lossy(&name),
            strand,
            String::from_utf8_lossy(&record.chrom),
            record.start,
            record.end
        );
    }
    Ok(())
}
```

### Reading a BED12 File

BED12 is the full format with exon information for gene structures:

```rust,no_run
use genepred::{Reader, Bed12};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = Reader::<Bed12>::from_path("data/transcripts.bed")?;
    for record in reader.records() {
        let record = record?;
        let name = record.name.unwrap_or(b"<unknown>".to_vec());
        println!(
            "Transcript: {} with {} exons",
            String::from_utf8_lossy(&name),
            record.exon_count()
        );
        for (i, (exon_start, exon_end)) in record.exons().iter().enumerate() {
            println!("  Exon {}: {}-{}", i + 1, exon_start, exon_end);
        }
    }
    Ok(())
}
```

### Using the Builder Pattern

For more control over reading behavior:

```rust,no_run
use genepred::{Reader, Bed6, ReaderMode};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = Reader::<Bed6>::builder()
        .from_path("data/large_file.bed")
        .mode(ReaderMode::Default)
        .buffer_capacity(128 * 1024)  // 128 KB buffer
        .build()?;

    for record in reader.records() {
        let record = record?;
        // Process record...
    }
    Ok(())
}
```

### Reading Files with Additional Custom Fields

If your BED file has extra columns beyond the standard format:

```rust,no_run
use genepred::{Reader, Bed6};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // BED6 + 2 additional fields
    let mut reader = Reader::<Bed6>::builder()
        .from_path("data/custom.bed")
        .additional_fields(2)
        .build()?;

    for record in reader.records() {
        let record = record?;
        let name = record.name.unwrap_or(b"<unknown>".to_vec());
        println!("Gene: {}", String::from_utf8_lossy(&name));

        // Additional fields are stored in extras with numeric keys
        if let Some(field1) = record.extras.get(&b"7".to_vec()) {
            println!("Custom field 1: {}", field1);
        }
        if let Some(field2) = record.extras.get(&b"8".to_vec()) {
            println!("Custom field 2: {}", field2);
        }
    }
    Ok(())
}
```

### Reading GTF/GFF Files

Parse GTF/GFF files into `GenePred` records. GTF/GFF readers work differently - 
they aggregate multiple rows into single GenePred records:

```rust,no_run
use genepred::{Reader, Gtf};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = Reader::<Gtf>::from_path("data/annotations.gtf")?;
    for record in reader.records() {
        let record = record?;
        println!(
            "Gene: {} with {} exons",
            String::from_utf8_lossy(&record.name.unwrap_or(b"<unknown>".to_vec())),
            record.exon_count()
        );
    }
    Ok(())
}
```

```rust,no_run
use genepred::{Reader, Gff};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = Reader::<Gff>::from_path("data/annotations.gff")?;
    for record in reader.records() {
        let record = record?;
        println!(
            "Gene: {} on {} strand",
            String::from_utf8_lossy(&record.name.unwrap_or(b"<unknown>".to_vec())),
            record.strand.map(|s| s.to_string()).unwrap_or_else(|| ".".to_string())
        );
    }
    Ok(())
}
```

**Note:** GTF and GFF formats now support the unified `from_path()` method alongside the traditional `from_gxf()` method. Both approaches provide the same functionality.

### Memory-Mapped Files for Maximum Performance

Use memory mapping for extremely fast access to large files:

```rust,no_run
use genepred::{Reader, Bed6, ReaderMode};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Enable the "mmap" feature in Cargo.toml
    let mut reader = Reader::<Bed6>::builder()
        .from_path("data/huge_file.bed")
        .mode(ReaderMode::Mmap)
        .build()?;

    for record in reader.records() {
        let record = record?;
        println!(
            "Record: {}",
            String::from_utf8_lossy(&record.name.unwrap_or(b"<unknown>".to_vec()))
        );
    }
    Ok(())
}
```

#### Direct Memory-Mapped Functions

For convenience, you can also use direct `from_mmap()` functions:

```rust,no_run
use genepred::{Reader, Bed6, Gtf, Gff};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // For BED files
    let mut bed_reader = Reader::<Bed6>::from_mmap("data/genes.bed")?;

    // For BED files with additional fields
    let mut bed_custom = Reader::<Bed6>::from_mmap_with_additional_fields("data/custom.bed", 2)?;

    // For GTF files
    let mut gtf_reader = Reader::<Gtf>::from_mmap("data/annotations.gtf")?;

    // For GFF files  
    let mut gff_reader = Reader::<Gff>::from_mmap("data/annotations.gff")?;

    // All work the same way
    for record in bed_reader.records() {
        let record = record?;
        println!("BED: {}", String::from_utf8_lossy(&record.name.unwrap_or(b"<unknown>".to_vec())));
    }

    for record in gtf_reader.records() {
        let record = record?;
        println!("GTF: {} exons", record.exon_count());
    }

    Ok(())
}
```

### Parallel Processing with Rayon

Process records in parallel for multi-core performance:

```rust,no_run
// Enable the "rayon" feature in Cargo.toml
use genepred::{Reader, Bed6};
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reader = Reader::<Bed6>::from_path("data/huge_file.bed")?;
    let counter = AtomicUsize::new(0);

    if let Ok(records) = reader.par_records() {
        records.for_each(|record| {
            match record {
                Ok(record) => {
                    // Process record in parallel
                    println!(
                        "Record: {}",
                        String::from_utf8_lossy(&record.name.unwrap_or(b"<unknown>".to_vec()))
                    );
                    counter.fetch_add(1, Ordering::Relaxed);
                }
                Err(e) => eprintln!("Error: {}", e),
            }
        });
    }
    println!("Processed {} records", counter.load(Ordering::Relaxed));
    Ok(())
}
```

#### Parallel Processing with Memory-Mapped Files

For maximum performance, combine memory mapping with parallel processing:

```rust,no_run
// Enable both "rayon" and "mmap" features in Cargo.toml
use genepred::{Reader, Bed6, Gtf};
use rayon::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parallel processing of BED files
    let bed_reader = Reader::<Bed6>::from_mmap("data/large.bed")?;
    if let Ok(records) = bed_reader.par_records() {
        let count = records
            .filter_map(Result::ok)
            .filter(|r| r.strand.map(|s| s.is_plus()).unwrap_or(false))
            .count();
        println!("Found {} records on plus strand", count);
    }

    // Parallel processing of GTF files
    let gtf_reader = Reader::<Gtf>::from_mmap("data/annotations.gtf")?;
    if let Ok(records) = gtf_reader.par_records() {
        let total_exons: usize = records
            .filter_map(Result::ok)
            .map(|r| r.exon_count())
            .sum();
        println!("Total exons: {}", total_exons);
    }

    Ok(())
}
```

### Reading Compressed Files

Automatically handle gzip-compressed files:

```rust,no_run
use genepred::{Reader, Bed3};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Enable the "compression" feature in Cargo.toml
    // Compression is auto-detected from .gz extension
    let mut reader = Reader::<Bed3>::from_path("data/regions.bed.gz")?;

    for record in reader.records() {
        let record = record?;
        println!("{:?}", record);
    }
    Ok(())
}
```

### Reading from Standard Input

Process streaming data from pipes:

```rust,no_run
use genepred::{Reader, Bed3};
use std::io;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let stdin = io::stdin();
    let mut reader = Reader::<Bed3>::from_reader(stdin)?;

    for record in reader.records() {
        let record = record?;
        // Process streaming records...
        println!("{:?}", record);
    }
    Ok(())
}
```

### Parallel I/O operations

I/O heavily parallelized

```rust,no_run
use genepred::{Bed12, Bed6, Gff, Reader, Writer};
use rayon::prelude::*;
use tempfile::NamedTempFile;

use std::collections::HashMap;
use std::path::PathBuf;

fn main() {
    let reader = Reader::<Bed12>::from_mmap(path).unwrap_or_else(|e| panic!("{}", e));
    let chunks = reader
        .par_chunks(10000)
        .unwrap_or_else(|e| panic!("{}", e))
        .map(|(idx, chunk)| {
            println!("Processing chunk {}", idx);
            let mut tmp = NamedTempFile::new().unwrap_or_else(|e| panic!("{}", e));

            chunk.into_iter().filter_map(Result::ok).for_each(|gene| {
                let _ = Writer::<Gff>::from_record(&gene, &mut tmp);
            });

            let (file, path) = tmp.keep().unwrap_or_else(|e| panic!("{}", e));
            path
        })
        .collect::<Vec<_>>();

    let mut writer = std::io::BufWriter::new(std::fs::File::create("output.gff").unwrap());

    for path in chunks {
        let mut file = std::fs::File::open(&path).unwrap_or_else(|e| panic!("{}", e));
        std::io::copy(&mut file, &mut writer).unwrap_or_else(|e| panic!("{}", e));

        std::fs::remove_file(path).unwrap_or_else(|e| panic!("{}", e));
    }
}
```


