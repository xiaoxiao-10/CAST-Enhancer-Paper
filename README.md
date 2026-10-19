# CAST-Enhancer — Paper Evidence Package

This is a **private manuscript-support repository** for the paper:

> **Nucleotide composition accounts for the shared component of enhancer activity across unseen cell types**

It is intentionally not a mirror of the old experimental workspace. The old
repository contains protocol history, exploratory analyses, raw-data staging,
failed model branches and intermediate checkpoints. This package keeps only the
files needed to audit and rebuild the current manuscript.

## Main result

The paper asks whether a pretrained genomic language model adds transferable
sequence information beyond nucleotide composition for FANTOM5 CAGE-defined
enhancer activity.

The answer is baseline-dependent:

- Against **logistic composition**, NT-LoRA looks better
  (+0.0175 AUPRC over all unseen loci; +0.0810 on unseen
  development-active loci).
- Against **XGBoost on the same 20 composition features**, NT-LoRA is
  indistinguishable on unseen development-active loci:
  **0.2415 versus 0.2425 AUPRC**, Δ = **−0.0011**
  (95% CI **−0.0358 to +0.0365**), NT better in **6/11** cell types.
- Adding 879 JASPAR motif counts does not improve on composition-XGBoost
  (0.2115 versus 0.2425).
- NT retains an advantage only on **retained** development-active loci
  (+0.0389; 95% CI 0.0155–0.0614), consistent with a locus-specific rather
  than transferable signal.
- No model improves on composition-XGBoost for development-inactive loci
  (NT − composition-XGBoost = −0.0020; 95% CI −0.0138 to +0.0007).

The manuscript therefore reports a **benchmark and baseline-dependence
finding**, not a new state-of-the-art enhancer model.

## Contents

| Directory | Purpose |
|---|---|
| `manuscript/en`, `manuscript/zh` | Final English and Chinese Word, PDF and Markdown manuscripts. |
| `figures/main` | Main Figures 1–4 in PNG, PDF and SVG. |
| `figures/supplementary` | Supplementary Figures S1–S3 in PNG, PDF and SVG. |
| `tables` | Main Tables 1–3 and Supplementary Tables S1–S8. |
| `results/locus_holdout` | Frozen per-cell, per-seed, macro, bootstrap, control and sensitivity tables. |
| `code/build` | Scripts that rebuild the English and Chinese Word manuscripts. |
| `code/figures` | Scripts that rebuild the main and supplementary figures. |
| `code/stats` | Archival post-processing and bootstrap scripts. |
| `code/training` | Archival training scripts, kept for protocol provenance. |
| `data/manifests` | Raw-data metadata, checksums and stage manifests. |
| `data/derived` | Derived inputs needed by the UpSet figure. |
| `docs` | Evidence map, reproducibility notes and data/code availability. |
| `provenance` | Source commit, environment and package checksums. |

## Reproduce the manuscript files

From the repository root:

```bash
python code/build/build_docx.py
python code/build/build_docx_zh.py
```

The Word files are written to `manuscript/en/` and `manuscript/zh/`.
The PDF files are exported from those Word files with Microsoft Word.

## Rebuild the figures

```bash
python code/figures/fig1_design_v10.py
python code/figures/fig2_baselines_v8.py
python code/figures/fig3_xgb_v8.py
python code/figures/fig4_control.py
python code/figures/supp_upset.py
python code/figures/figS_overlap.py
```

The scripts expect the frozen aggregate tables in `results/locus_holdout/`.
Training scripts under `code/training/` contain the original remote-path
assumptions and are archival; they are not required to rebuild the manuscript
or the figures.

## What is not included

- Raw FANTOM5, ENCODE, hg38 or JASPAR files.
- Model checkpoints, prediction arrays or zarr stores.
- Temporary render directories, failed branches or exploratory notebooks.
- Any data that cannot be redistributed under the upstream licence.

## Data and code availability

Upstream data are public: the FANTOM5 hg38 enhancer atlas, ENCODE K562,
HepG2 and GM12878 ATAC-seq, the hg38 primary assembly and JASPAR 2024.

This private package includes aggregate result tables and the code used to
rebuild the manuscript figures and documents. **The repository is not a public
code or weight release.** Model weights and the full training pipeline remain
available from the corresponding author on reasonable request for
non-commercial research, subject to the base model licence.

See `docs/DATA_AVAILABILITY.md` for the exact wording used in the manuscript.

## Provenance

The package was migrated from the `paper/v3-manuscript-draft` branch of the
old CAST-Enhancer workspace at commit `7fb20a2`.
See `provenance/SOURCE_COMMIT.txt` and `provenance/SHA256SUMS.txt`.
