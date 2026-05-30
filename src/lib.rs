//! FPKM computation from BAM + BED12 gene model.
//!
//! Formula: FPKM = (Frag_count × 10⁹) / (mRNA_size × total_fragments)
//!          FPM  = (Frag_count / total_fragments) × 10⁶
//!
//! where mRNA_size = sum of exon block lengths from BED12 columns 10/11,
//! and total_fragments is the count of SE reads plus read1 of PE pairs.
//!
//! Fragment-to-gene assignment follows RSeQC's counting convention:
//! - Single-end reads: the full read span is checked against each exon interval.
//! - Paired-end reads: only read1 is processed; read2 is skipped entirely.
//!   Overlap is tested as a point check on the read1 start position (1 bp), not
//!   the full span — this matches FPKM_count.py's `exon_ranges.find(pos, pos+1)`.
//!
//! Reads/pairs hitting 2+ distinct genes are discarded as ambiguous.

#![allow(clippy::cast_precision_loss)]

pub mod bed;
pub mod count;
pub mod output;
pub mod strand;

pub use bed::{Gene, load_bed12};
pub use count::{CountResult, FpkmOpts, count_bam};
pub use output::write_fpkm;
pub use strand::StrandRules;
