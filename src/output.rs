use std::io::{BufWriter, Write};

use rsomics_common::{Result, RsomicsError};

use crate::bed::Gene;
use crate::count::CountResult;

/// Write the FPKM table to `output`.
///
/// Tab-separated with a `#`-prefixed header line:
/// `#chrom\tst\tend\taccession\tmRNA_size\tgene_strand\tFrag_count\tFPM\tFPKM`
pub fn write_fpkm<W: Write>(genes: &[Gene], result: &CountResult, output: W) -> Result<()> {
    let mut w = BufWriter::with_capacity(256 * 1024, output);
    writeln!(
        w,
        "#chrom\tst\tend\taccession\tmRNA_size\tgene_strand\tFrag_count\tFPM\tFPKM"
    )
    .map_err(RsomicsError::Io)?;

    let total = result.total_fragments;
    for (gene, &count) in genes.iter().zip(result.frag_counts.iter()) {
        let fpm = if total > 0.0 {
            count / total * 1_000_000.0
        } else {
            0.0
        };
        let fpkm = if total > 0.0 && gene.mrna_size > 0 {
            count * 1_000_000_000.0 / (gene.mrna_size as f64 * total)
        } else {
            0.0
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            gene.chrom,
            gene.tx_start,
            gene.tx_end,
            gene.name,
            py_float(gene.mrna_size as f64),
            gene.strand,
            py_float(count),
            py_float(fpm),
            py_float(fpkm),
        )
        .map_err(RsomicsError::Io)?;
    }

    w.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

/// Python-compatible float rendering: integer values get the `.0` suffix.
///
/// RSeQC uses Python's `str(float)` which renders `2.0` not `2`.
fn py_float(v: f64) -> String {
    if v.fract() == 0.0 && v.is_finite() {
        format!("{v}.0")
    } else {
        format!("{v}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bed::Gene;
    use crate::count::CountResult;

    #[test]
    fn fpkm_formula() {
        // Frag_count=2, total=4, mrna_size=250
        // FPM = 2/4 * 1e6 = 500000; FPKM = 2e9/(250*4) = 2000000
        let genes = vec![Gene {
            chrom: "chr1".into(),
            tx_start: 100,
            tx_end: 500,
            name: "G1".into(),
            strand: '+',
            mrna_size: 250,
            exons: vec![(100, 200), (350, 500)],
        }];
        let result = CountResult {
            frag_counts: vec![2.0],
            total_fragments: 4.0,
        };
        let mut buf = Vec::new();
        write_fpkm(&genes, &result, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("500000"), "{s}");
        assert!(s.contains("2000000"), "{s}");
    }
}
