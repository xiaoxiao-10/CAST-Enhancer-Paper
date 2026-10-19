#!/usr/bin/env bash
set -euo pipefail

OUT="/home/vipuser/cast-v2-locus-holdout"
PY="/home/vipuser/miniforge3/envs/cast-enhancer/bin/python"

while [ "$(find "$OUT" -maxdepth 1 -name 'LH_scores_ntinit_seed_*.npy' | wc -l)" -lt 5 ]; do
  sleep 60
done

cd /home/vipuser/cast-v2-v0
"$PY" lh_ntinit_posthoc.py \
  --out-dir "$OUT" \
  --require-seeds 5 \
  --boot 1000 \
  --macro-boot 3000 \
  > "$OUT/ntinit_posthoc.log" 2>&1

touch "$OUT/ntinit_posthoc.done"
