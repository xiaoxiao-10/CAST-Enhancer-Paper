# S2. Protocol history and correction record

## 1. Full-locus model (2026-09-26)

The original v2 analysis trained on all quality-passing loci and evaluated the
same loci in 11 held-out primary cell types. The five-seed ensemble reached a
macro lift of 10.333× for NT-LoRA, 6.017× for the composition baseline and
5.419× for the DNA-CNN (`v2/UPDATE_20260926_FINAL.md`,
`v2/V15_final_heldout_analysis.tsv`).

Because training and evaluation shared the same candidate loci, this design
could not separate sequence generalisation from locus-level recognition.

## 2. Chromosome-level audit (2026-10-17)

An audit split the full-locus model's held-out-cell results by chromosome
provenance. On development-active loci, the NT-minus-composition gain was
+0.084 on training chromosomes and +0.006 on held-out chromosomes, and 54 of 55
cell-by-seed pairs showed a larger gain on training chromosomes
(`v2/audit_20261017/`).

That analysis had two limitations. The held-out-chromosome subset contained
about 55 positives per cell, giving a hierarchical 95% interval of −0.032 to
+0.044; and the direction reversed on unseen loci in the subsequent experiment.

## 3. Locus holdout (2026-10-18)

To test locus memorisation directly, 20% of quality-passing loci (12,360 of
61,752, stratified by chromosome and GC) were excluded from all development-cell
training and all three model families were retrained on the remaining 80%
(108,252 training, 12,621 validation and 24,247 test records).

Five-seed results:

| Locus set | Development activity | Composition | CNN | NT | NT − composition (95% CI) |
|---|---|---:|---:|---:|---|
| unseen | all | 5.66× | 5.06× | 8.35× | +0.0175 (+0.0101, +0.0291) |
| unseen | active | 3.96× | 3.69× | 5.94× | +0.0810 (+0.0501, +0.1219) |
| retained | active | 5.70× | 4.96× | 6.77× | +0.0487 (+0.0203, +0.0747) |
| unseen | inactive | 2.74× | 2.55× | 2.22× | −0.0007 (−0.0073, +0.0013) |

The point estimate was positive in all 11 held-out cell types on unseen
development-active loci (range +0.0280 to +0.1719); per-cell 95% intervals
excluded zero in 6 of 11. The point estimate was larger on unseen than on
retained loci, but the difference was not resolved (+0.0323; 95% CI −0.0036 to
+0.0768; bootstrap p = 0.074). The values in the table are point estimates;
`LH_posthoc_macro.tsv` retains the bootstrap means used to construct intervals.

## 4. Correction

The chromosome-level reading in section 2 — that the advantage was confined to
training chromosomes — is not supported by the locus holdout. It is retained as
a sensitivity analysis (`Fig. 4`) and is not used as a primary claim. The
primary generalisation test is the locus holdout.

## 5. Random-init control (2026-10-18)

The same NT architecture was initialised from the model config without
pretrained weights and trained with the same tokenizer, LoRA modules, data,
optimizer and four-epoch schedule for the same five seeds. On retained
development-active loci, the random-init ensemble reached AUPRC 0.258 versus
0.273 for the pretrained ensemble (Δ +0.0153; 95% CI −0.0147 to +0.0465). On
unseen development-active loci, it fell to 0.114 versus 0.241 (Δ +0.1275; 95%
CI +0.0917 to +0.1812; positive in 11/11 cell types), and was below composition
(−0.0465; 95% CI −0.0723 to −0.0303). The retained-to-unseen dissociation is
consistent with locus-level memorization in the random-init model and a
pretraining-dependent generalization advantage. The control uses the same
four-epoch budget and was not separately tuned.

## 6. Methodological conclusion

Holding out cell types alone does not establish sequence generalisation when
candidate loci are shared between training and evaluation. A locus or chromosome
holdout is required in addition. In this study, the NT gain transfers to unseen
loci, but only within the shared activity component; loci that are inactive in
every development cell are not predicted better than by composition. The
random-init control shows that the unseen-locus advantage depends on pretrained
weights under the tested training budget.
