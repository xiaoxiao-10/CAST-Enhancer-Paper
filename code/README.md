# Code notes

The code is split by function:

- `code/build/` rebuilds the English and Chinese Word manuscripts.
- `code/figures/` rebuilds the main and supplementary figures.
- `code/stats/` contains the post-processing and bootstrap scripts.
- `code/training/` contains the archival training scripts.

The figure and manuscript-build scripts have been adapted to this package's
layout and are runnable from the repository root. The training and some
post-processing scripts retain the original remote Linux paths
(`/home/vipuser/...`), because they document the exact analysis run rather
than a portable training pipeline.

The frozen aggregate tables in `results/locus_holdout/` are sufficient to
rebuild every manuscript number and every figure.
