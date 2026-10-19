# Data and code availability

## Manuscript wording

> Code and trained weights are not publicly released. They are available from
> the corresponding author on reasonable request for non-commercial research,
> subject to the licence of the base Nucleotide Transformer model. Aggregate
> outputs needed to reproduce the reported tables and figures, including cohort
> counts, per-cell metrics, per-seed metrics, the CD19-exclusion sensitivity
> analysis and the sequence-neighborhood audit, are provided as Supplementary
> Tables S1–S8. Upstream data are the FANTOM5 hg38 enhancer atlas, the hg38
> primary assembly and the ENCODE hg38 blacklist v2, with per-file URLs and
> checksums in the supplementary data package.

## Upstream resources

- FANTOM5 hg38 enhancer atlas:
  `F5.hg38.enhancers.bed.gz`,
  `F5.hg38.enhancers.expression.usage.matrix.gz`,
  `F5.hg38.enhancers.expression.tpm.matrix.gz`.
- ENCODE K562, HepG2 and GM12878 ATAC-seq.
- UCSC hg38 primary assembly.
- ENCODE hg38 blacklist v2 (`ENCFF356LFX`).
- JASPAR 2024 CORE vertebrates non-redundant PFMs.

## Repository boundary

This is a private manuscript-support package. It contains aggregate result
tables and manuscript/figures code, but not raw data or model weights. Do not
make the repository public without author approval.
