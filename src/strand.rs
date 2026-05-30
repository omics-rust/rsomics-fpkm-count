use rsomics_common::{Result, RsomicsError};

/// Strand specificity rules parsed from RSeQC's `-d` format.
///
/// Each token is two or three characters:
/// - Optional leading `1` (read1) or `2` (read2); absent → applies to both.
/// - Read strand: `+` = forward, `-` = reverse.
/// - Gene strand: `+` = forward, `-` = reverse.
///
/// Examples: `"1++,1--,2+-,2-+"`, `"++,--"`, `"+-,-+"`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct StrandRules {
    /// Bitmask: bit index = (is_read2 as u8)*4 + (read_reverse as u8)*2 + (gene_minus as u8).
    mask: u8,
}

impl StrandRules {
    const fn bit(is_read2: bool, read_reverse: bool, gene_minus: bool) -> u8 {
        1 << ((is_read2 as u8) * 4 + (read_reverse as u8) * 2 + gene_minus as u8)
    }

    pub fn allows(&self, is_read2: bool, read_reverse: bool, gene_strand: char) -> bool {
        self.mask & Self::bit(is_read2, read_reverse, gene_strand == '-') != 0
    }

    pub fn parse(s: &str) -> Result<Option<Self>> {
        let s = s.trim();
        if s.is_empty() {
            return Ok(None);
        }
        let mut rules = Self::default();
        for token in s.split(',') {
            let token = token.trim();
            let (r2, rest) = if let Some(t) = token.strip_prefix('1') {
                (Some(false), t)
            } else if let Some(t) = token.strip_prefix('2') {
                (Some(true), t)
            } else {
                (None, token)
            };
            if rest.len() != 2 {
                return Err(RsomicsError::InvalidInput(format!(
                    "strand rule token {token:?}: expected 2-char strand pair after optional read number"
                )));
            }
            let read_fwd = rest.as_bytes()[0] == b'+';
            let gene_fwd = rest.as_bytes()[1] == b'+';
            match r2 {
                Some(is_r2) => {
                    rules.mask |= Self::bit(is_r2, !read_fwd, !gene_fwd);
                }
                None => {
                    rules.mask |= Self::bit(false, !read_fwd, !gene_fwd);
                    rules.mask |= Self::bit(true, !read_fwd, !gene_fwd);
                }
            }
        }
        Ok(Some(rules))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strand_rules_parse_paired_end() {
        let r = StrandRules::parse("1++,1--,2+-,2-+").unwrap().unwrap();
        assert!(r.allows(false, false, '+'));
        assert!(r.allows(false, true, '-'));
        assert!(r.allows(true, false, '-'));
        assert!(r.allows(true, true, '+'));
        assert!(!r.allows(false, false, '-'));
        assert!(!r.allows(true, false, '+'));
    }

    #[test]
    fn strand_rules_parse_single_end() {
        let r = StrandRules::parse("++,--").unwrap().unwrap();
        assert!(r.allows(false, false, '+'));
        assert!(r.allows(false, true, '-'));
        assert!(!r.allows(false, false, '-'));
    }

    #[test]
    fn strand_rules_none_on_empty() {
        assert!(StrandRules::parse("").unwrap().is_none());
    }
}
