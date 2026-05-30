use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use noodles::bam;
use rsomics_common::{Result, RsomicsError};

use crate::bed::Gene;
use crate::strand::StrandRules;

pub struct FpkmOpts {
    /// Minimum MAPQ to be considered uniquely mapped.
    pub min_mapq: u8,
    /// Skip reads with MAPQ below `min_mapq` (mirrors `-u`).
    pub unique_only: bool,
    /// Strand specificity rules; `None` = unstranded.
    pub strand: Option<StrandRules>,
}

impl Default for FpkmOpts {
    fn default() -> Self {
        Self {
            min_mapq: 30,
            unique_only: false,
            strand: None,
        }
    }
}

/// Fragment count per gene and total aligned fragments.
pub struct CountResult {
    /// Fragment count per gene, parallel to the input `genes` slice.
    pub frag_counts: Vec<f64>,
    /// Total unique aligned fragments (pairs deduplicated by QNAME).
    pub total_fragments: f64,
}

/// Count fragments per gene from `bam_path`.
pub fn count_bam(bam_path: &Path, genes: &[Gene], opts: &FpkmOpts) -> Result<CountResult> {
    let exon_by_chrom = build_exon_index(genes);

    let file = File::open(bam_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bam_path.display())))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().map_err(RsomicsError::Io)?;

    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(ToString::to_string)
        .collect();

    let mut frag_counts = vec![0_f64; genes.len()];
    let mut total_fragments: f64 = 0.0;

    for result in reader.records() {
        let record = result.map_err(RsomicsError::Io)?;
        let flags = record.flags();

        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        let mq = record.mapping_quality().map_or(0, |q| q.get());
        if opts.unique_only && mq < opts.min_mapq {
            continue;
        }

        let is_paired = flags.is_segmented();
        let is_read2 = flags.is_last_segment();

        // RSeQC counts total_fragments as: SE reads + read1 of PE pairs.
        // read2 is skipped entirely for both total counting and gene assignment.
        if is_read2 {
            continue;
        }

        total_fragments += 1.0;

        let Some(tid) = record.reference_sequence_id().transpose().ok().flatten() else {
            continue;
        };
        let Some(chrom) = ref_names.get(tid) else {
            continue;
        };
        let Some(pos) = record.alignment_start().transpose().ok().flatten() else {
            continue;
        };
        let read_start = pos.get() as u64 - 1;
        // SE: check full read span; PE read1: point check on read start (RSeQC convention).
        let overlap_end = if is_paired {
            read_start + 1
        } else {
            read_start + record.sequence().len() as u64
        };
        let read_reverse = flags.is_reverse_complemented();

        let hits = exon_gene_hits(
            chrom,
            read_start,
            overlap_end,
            is_read2,
            read_reverse,
            genes,
            &exon_by_chrom,
            opts,
        );

        if hits.len() == 1 {
            frag_counts[hits[0]] += 1.0;
        }
    }

    Ok(CountResult {
        frag_counts,
        total_fragments,
    })
}

fn build_exon_index(genes: &[Gene]) -> HashMap<String, Vec<(u64, u64, usize)>> {
    let mut m: HashMap<String, Vec<(u64, u64, usize)>> = HashMap::new();
    for (gi, gene) in genes.iter().enumerate() {
        for &(es, ee) in &gene.exons {
            m.entry(gene.chrom.clone()).or_default().push((es, ee, gi));
        }
    }
    m
}

#[allow(clippy::too_many_arguments)]
fn exon_gene_hits(
    chrom: &str,
    read_start: u64,
    read_end: u64,
    is_read2: bool,
    read_reverse: bool,
    genes: &[Gene],
    exon_by_chrom: &HashMap<String, Vec<(u64, u64, usize)>>,
    opts: &FpkmOpts,
) -> Vec<usize> {
    let Some(exons) = exon_by_chrom.get(chrom) else {
        return Vec::new();
    };
    let mut hits: Vec<usize> = Vec::new();
    for &(es, ee, gi) in exons {
        if read_start >= ee || read_end <= es {
            continue;
        }
        if hits.contains(&gi) {
            continue;
        }
        if let Some(rules) = opts.strand
            && !rules.allows(is_read2, read_reverse, genes[gi].strand)
        {
            continue;
        }
        hits.push(gi);
    }
    hits
}
