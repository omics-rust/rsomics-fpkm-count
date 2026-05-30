//! FPKM/FPM computation from BAM + BED12.
//!
//! FPKM = (count × 10⁹) / (mRNA_size × total_fragments); mRNA_size = sum of BED12 exon blocks.
//! Assignment follows RSeQC: SE uses the full span; PE read1 uses a 1 bp point check
//! at read start (FPKM_count.py `exon_ranges.find(pos, pos+1)`). Reads hitting
//! 2+ genes are discarded as ambiguous.

#![allow(clippy::cast_precision_loss)]

pub mod bed;
pub mod count;
pub mod output;
pub mod strand;

pub use bed::{Gene, load_bed12};
pub use count::{CountResult, FpkmOpts, count_bam};
pub use output::write_fpkm;
pub use strand::StrandRules;
