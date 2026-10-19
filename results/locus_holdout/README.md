# Locus-holdout results

Frozen aggregate outputs for the paper. The manuscript uses the
`LH_correct_*` tables as the current label-eligible analysis.

## Core files

| File | Content |
|---|---|
| `LH_correct_macro.tsv` | Macro AUPRC for logistic composition, composition-XGBoost, motif-XGBoost, composition+motif-XGBoost and NT-LoRA in six subset/activity strata. |
| `LH_correct_boot.tsv` | Macro 95% paired 1 Mb block-bootstrap intervals for NT − composition-XGBoost and NT − composition+motif-XGBoost. |
| `LH_correct_metrics.tsv` | Per-cell AUPRC for the five models and the two stronger-baseline deltas. |
| `LH_model_compare_macro.tsv` | Long-form macro AUPRC for the stronger baselines. |
| `LH_model_compare_metrics.tsv` | Long-form per-cell AUPRC for the stronger baselines. |
| `LH_nt_vs_strong_macro.tsv` | NT versus composition-XGBoost and composition+motif-XGBoost, macro. |
| `LH_nt_vs_strong_boot.tsv` | Bootstrap intervals for the stronger-baseline deltas. |
| `LH_nt_vs_strong_baselines.tsv` | Per-cell NT and baseline AUPRC with deltas. |

## Controls and sensitivity

| File | Content |
|---|---|
| `LH_ntinit_macro.tsv` | Random-init NT versus pretrained NT and composition, macro. |
| `LH_ntinit_deltas.tsv` | Per-cell random-init deltas with intervals. |
| `LH_ntinit_metrics.tsv` | Per-cell AUPRC for composition, pretrained NT and random-init NT. |
| `LH_nocell_macro.tsv` | Cell-free NT control, macro. |
| `LH_nocell_deltas.tsv` | Cell-free NT control, per-cell deltas. |
| `LH_nocell_metrics.tsv` | Cell-free NT control, per-cell metrics. |
| `LH_overlap_macro.tsv` | Sequence-disjoint overlap sensitivity, macro. |
| `LH_overlap_deltas.tsv` | Sequence-disjoint overlap sensitivity, per-cell deltas. |
| `LH_overlap_metrics.tsv` | Sequence-disjoint overlap sensitivity, per-cell metrics. |
| `LH_chrom_sensitivity.tsv` | Chromosome-provenance sensitivity. |
| `LH_posthoc_*` | Earlier logistic-composition analysis, retained for provenance and the NT-versus-logistic comparison. |

## Current conclusion

On unseen development-active loci:

- composition-XGBoost: **0.2425**
- NT-LoRA: **0.2415**
- NT − composition-XGBoost: **−0.0011**
  (95% CI **−0.0358 to +0.0365**)
- NT better in **6/11** held-out cell types

On retained development-active loci, NT − composition-XGBoost is
**+0.0389** (95% CI **0.0155 to 0.0614**). This is consistent with a
locus-specific rather than transferable signal.

On unseen development-inactive loci, NT − composition-XGBoost is
**−0.0020** (95% CI **−0.0138 to +0.0007**). No model improves on the
strongest composition baseline for the cell-specific component.

## Note on table versions

- `LH_correct_*` is the current analysis.
- `LH_posthoc_*` is the earlier logistic-composition analysis.
- Do not use old 0-filled UpSet tables for AUPRC or prevalence.
