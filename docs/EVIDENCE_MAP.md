# Evidence map

Every main-text claim should map to a frozen table or figure in this package.
The manuscript is written against the **label-eligible** locus-holdout tables,
not against the earlier 0-filled UpSet tables.

| Claim | Manuscript location | Primary evidence |
|---|---|---|
| 63,285 FANTOM5 candidate loci; 20% locus holdout; 11 held-out primary cell types | Results 2.1; Table 1 | `tables/table1_cohort.tsv`; `results/locus_holdout/LH_dev_split_counts.tsv`; `results/locus_holdout/heldout_label_counts.tsv` |
| 21-parameter logistic composition baseline reaches 5.66× on unseen loci | Results 2.2; Table 2 | `tables/table2_four_strata.tsv`; `results/locus_holdout/LH_correct_macro.tsv` |
| DNA-CNN does not exceed the composition baseline | Results 2.3; Table 2 | `tables/table2_four_strata.tsv`; `results/locus_holdout/LH_posthoc_metrics.tsv` |
| NT exceeds logistic composition on unseen loci (+0.0175; +0.0810 within development-active loci) | Results 2.4; Table 2; Fig. 3 (legacy context) | `results/locus_holdout/LH_macro_point.tsv`; `results/locus_holdout/LH_posthoc_macro.tsv`; `results/locus_holdout/LH_posthoc_deltas.tsv` |
| Composition-XGBoost reaches 0.243 on unseen development-active loci, versus 0.241 for NT | Results 2.6; Table 3 | `tables/table3_strong_baselines.tsv`; `results/locus_holdout/LH_correct_macro.tsv`; `results/locus_holdout/LH_model_compare_macro.tsv` |
| NT − composition-XGBoost = −0.0011 (95% CI −0.0358 to +0.0365); NT better in 6/11 cells | Results 2.6; Table 3; Fig. 3 | `results/locus_holdout/LH_correct_macro.tsv`; `results/locus_holdout/LH_correct_boot.tsv`; `results/locus_holdout/LH_correct_metrics.tsv` |
| Motif-only and composition+motif do not improve on composition-XGBoost | Results 2.6; Table 3; Discussion | `tables/table3_strong_baselines.tsv`; `results/locus_holdout/LH_model_compare_macro.tsv`; `results/locus_holdout/LH_motif_macro.tsv` |
| NT retains +0.0389 on retained development-active loci but not on unseen loci | Results 2.6; Table 3 | `tables/table3_strong_baselines.tsv`; `results/locus_holdout/LH_correct_macro.tsv`; `results/locus_holdout/LH_correct_boot.tsv` |
| No advantage is detected in the cell-specific (development-inactive) component | Results 2.7; Table 3 | `tables/table3_strong_baselines.tsv`; `results/locus_holdout/LH_correct_macro.tsv` |
| Chromosome-level splits are low-powered and unstable | Results 2.8; Supplementary Fig. S1 | `results/locus_holdout/LH_chrom_sensitivity.tsv`; `figures/supplementary/figS_chromosome.*` |
| Random-init NT matches pretrained NT on retained loci but falls below composition on unseen loci | Results 2.9; Fig. 4; Table S4 | `results/locus_holdout/LH_ntinit_macro.tsv`; `LH_ntinit_deltas.tsv`; `LH_ntinit_metrics.tsv`; `tables/supplementary/S4_ntinit_macro.tsv` |
| Cell-free NT control still beats logistic composition but not composition-XGBoost | Results 2.9; Table S8 | `results/locus_holdout/LH_nocell_macro.tsv`; `tables/supplementary/S8_nocell_macro.tsv` |
| Sequence-disjoint sensitivity attenuates the logistic-baseline gain but does not create a composition-XGBoost advantage | Results 2.10; Table S6; Supplementary Fig. S3 | `results/locus_holdout/LH_overlap_macro.tsv`; `LH_overlap_deltas.tsv`; `LH_overlap_metrics.tsv`; `tables/supplementary/S6_overlap_sensitivity.tsv` |
| CD19+ B-cell exclusion does not change the conclusion | Limitations; Table S5 | `tables/supplementary/S5_no_cd19.tsv` |

## Important table distinction

- `LH_correct_*.tsv` is the current label-eligible analysis used in the
  manuscript.
- `LH_posthoc_*.tsv` is the earlier logistic-composition analysis, retained
  for the NT-versus-logistic comparison and provenance.
- Do not use the old 0-filled UpSet tables for prevalence or AUPRC claims.
- `LH_model_compare_macro.tsv` and `LH_nt_vs_strong_macro.tsv` are the source
  for the stronger-baseline Table 3.
