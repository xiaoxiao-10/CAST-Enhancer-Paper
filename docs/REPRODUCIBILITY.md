# Reproducibility notes

## Environment

The manuscript-building and figure scripts were run with:

- Python 3.11
- matplotlib 3.8
- pandas 2.x
- numpy 1.26
- python-docx 1.x
- Pillow 10.x
- PyMuPDF 1.24

See `requirements.txt` for a minimal environment.

## Rebuild the Word manuscripts

```bash
python code/build/build_docx.py
python code/build/build_docx_zh.py
```

The scripts read:

- `manuscript/en/manuscript_v3.md`
- `manuscript/zh/manuscript_v3_zh.md`
- `figures/main/*.png`
- `tables/table1_cohort.tsv`
- `tables/table2_four_strata.tsv`
- `results/locus_holdout/*.tsv`

and write the Word files into `manuscript/en/` and `manuscript/zh/`.

PDF export is a separate Microsoft Word step:

1. open the `.docx` file in Word;
2. export as PDF with bookmarks and document properties enabled;
3. keep the PDF names used in the package.

## Rebuild the figures

```bash
python code/figures/fig1_design_v10.py
python code/figures/fig2_baselines_v8.py
python code/figures/fig3_xgb_v8.py
python code/figures/fig4_control.py
python code/figures/supp_upset.py
python code/figures/figS_overlap.py
```

The scripts write PNG/PDF/SVG/TIFF into `figures/main` or
`figures/supplementary`. TIFF files are intentionally excluded from the
repository because they are large; they can be regenerated with the same
scripts.

## Training and post-processing

`code/training/` contains the original training scripts for provenance.
They assume the remote Linux layout (`/home/vipuser/...`) and are not
self-contained. `code/stats/` contains the post-processing and bootstrap
scripts; some of them also retain the original remote-path assumptions.

The frozen aggregate tables in `results/locus_holdout/` are sufficient to
rebuild every main-text number and every figure. Retraining the models is not
required to audit the manuscript.

## Reproducibility boundary

- This package is not a public code or weight release.
- Raw data and model checkpoints are not included.
- The frozen tables are the authoritative source for the manuscript numbers.
- Any change to labels, splits, filters or model outputs requires a new
  evidence version and a new manuscript number.
