use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// One BED12 transcript record with exon intervals pre-computed from block geometry.
#[derive(Debug, Clone)]
pub struct Gene {
    pub chrom: String,
    pub tx_start: u64,
    pub tx_end: u64,
    pub name: String,
    pub strand: char,
    /// Sum of BED12 block sizes — the effective transcript length for the FPKM denominator.
    pub mrna_size: u64,
    /// Exon intervals in BED half-open [start, end) coordinates.
    pub exons: Vec<(u64, u64)>,
}

/// Load a BED12 file and return one [`Gene`] per transcript record.
pub fn load_bed12(path: &Path) -> Result<Vec<Gene>> {
    let f = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut genes = Vec::new();
    for (i, line) in BufReader::new(f).lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        let s = line.trim_end();
        if s.is_empty() || s.starts_with('#') || s.starts_with("track") || s.starts_with("browser")
        {
            continue;
        }
        genes.push(parse_bed12_line(s, i + 1)?);
    }
    Ok(genes)
}

fn parse_bed12_line(line: &str, lineno: usize) -> Result<Gene> {
    let mut f = line.split('\t');
    macro_rules! field {
        ($name:expr) => {
            f.next().ok_or_else(|| {
                RsomicsError::InvalidInput(format!("BED12 line {lineno}: missing {}", $name))
            })?
        };
    }
    let chrom = field!("chrom").to_string();
    let tx_start: u64 = field!("chromStart")
        .parse()
        .map_err(|e| RsomicsError::InvalidInput(format!("BED12 line {lineno}: chromStart: {e}")))?;
    let tx_end: u64 = field!("chromEnd")
        .parse()
        .map_err(|e| RsomicsError::InvalidInput(format!("BED12 line {lineno}: chromEnd: {e}")))?;
    let name = field!("name").to_string();
    let _score = field!("score");
    let strand = field!("strand").chars().next().unwrap_or('.');
    let _thick_start = field!("thickStart");
    let _thick_end = field!("thickEnd");
    let _item_rgb = field!("itemRgb");
    let block_count: usize = field!("blockCount")
        .parse()
        .map_err(|e| RsomicsError::InvalidInput(format!("BED12 line {lineno}: blockCount: {e}")))?;
    let sizes_raw = field!("blockSizes");
    let starts_raw = field!("blockStarts");

    let sizes = comma_u64(sizes_raw, block_count, lineno, "blockSizes")?;
    let starts = comma_u64(starts_raw, block_count, lineno, "blockStarts")?;

    let mut exons = Vec::with_capacity(block_count);
    let mut mrna_size = 0u64;
    for i in 0..block_count {
        let es = tx_start + starts[i];
        let ee = es + sizes[i];
        mrna_size += sizes[i];
        exons.push((es, ee));
    }

    Ok(Gene {
        chrom,
        tx_start,
        tx_end,
        name,
        strand,
        mrna_size,
        exons,
    })
}

fn comma_u64(s: &str, expected: usize, lineno: usize, field: &str) -> Result<Vec<u64>> {
    let mut v = Vec::with_capacity(expected);
    for p in s.split(',') {
        let p = p.trim();
        if p.is_empty() {
            continue;
        }
        v.push(p.parse::<u64>().map_err(|e| {
            RsomicsError::InvalidInput(format!("BED12 line {lineno}: {field} value {p:?}: {e}"))
        })?);
    }
    if v.len() != expected {
        return Err(RsomicsError::InvalidInput(format!(
            "BED12 line {lineno}: {field} has {} values, expected {expected}",
            v.len()
        )));
    }
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bed12_parse_simple() {
        let line = "chr1\t100\t500\tGENE1\t0\t+\t100\t500\t0\t2\t100,150,\t0,250,";
        let gene = parse_bed12_line(line, 1).unwrap();
        assert_eq!(gene.name, "GENE1");
        assert_eq!(gene.tx_start, 100);
        assert_eq!(gene.tx_end, 500);
        assert_eq!(gene.mrna_size, 250);
        assert_eq!(gene.exons, vec![(100, 200), (350, 500)]);
    }
}
