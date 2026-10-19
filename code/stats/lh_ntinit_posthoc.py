#!/usr/bin/env python
"""Evaluate the random-initialised Nucleotide Transformer control.

The comparison is paired within cell type and locus stratum:
  nt      = pretrained NT-LoRA, five-seed ensemble
  ntinit  = same architecture and tokenizer, random initialisation, five-seed ensemble
  comp    = 21-parameter mono/dinucleotide logistic regression

Outputs are written next to the frozen locus-holdout results.  Existing files
are never overwritten; the script writes new LH_ntinit_* files.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

sys.path.insert(0, "/home/vipuser/cast-v2-v0")
import lh_train as L


SEEDS = [13, 37, 71, 101, 137]
KEY_STRATA = [
    ("unseen", "all"),
    ("unseen", "dev_active"),
    ("seen", "dev_active"),
    ("unseen", "dev_inactive"),
]


def load_seed_ensemble(out_dir, tag, require_seeds):
    files = [os.path.join(out_dir, "LH_scores_%s_seed_%d.npy" % (tag, s)) for s in SEEDS]
    present = [(s, f) for s, f in zip(SEEDS, files) if os.path.exists(f)]
    if len(present) < require_seeds:
        raise SystemExit(
            "only %d/%d %s seeds present: %s"
            % (len(present), require_seeds, tag, [os.path.basename(f) for _, f in present])
        )
    arrays = [np.load(f) for _, f in present]
    n = len(arrays[0])
    if any(len(a) != n for a in arrays):
        raise SystemExit("%s seed score files have inconsistent lengths" % tag)
    return np.mean(arrays, axis=0), [s for s, _ in present]


def block_indices(chrom, mid):
    blocks = np.array(
        [c + ":" + str(p) for c, p in zip(chrom, mid // 1000000)],
        dtype=object,
    )
    unique = np.unique(blocks)
    groups = [np.where(blocks == b)[0] for b in unique]
    return groups


def paired_block_delta(y, score_a, score_b, groups, rng, n_boot):
    """Return point delta, percentile CI, skipped count, and bootstrap draws."""
    point = average_precision_score(y, score_a) - average_precision_score(y, score_b)
    draws = []
    skipped = 0
    for _ in range(n_boot):
        picked = rng.integers(0, len(groups), len(groups))
        idx = np.concatenate([groups[i] for i in picked])
        if y[idx].sum() < 1:
            skipped += 1
            continue
        draws.append(
            average_precision_score(y[idx], score_a[idx])
            - average_precision_score(y[idx], score_b[idx])
        )
    if not draws:
        return point, np.nan, np.nan, skipped, np.array([], dtype=float)
    draws = np.asarray(draws, dtype=float)
    return (point, float(np.percentile(draws, 2.5)),
            float(np.percentile(draws, 97.5)), skipped, draws)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/home/vipuser/cast-v2-locus-holdout")
    ap.add_argument("--cnn-ensemble", default=None)
    ap.add_argument("--require-seeds", type=int, default=5)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--macro-boot", type=int, default=3000)
    args = ap.parse_args()

    out_dir = args.out_dir
    idx, _, _, _ = L.load_inputs()
    holdout = np.load(os.path.join(out_dir, "LH_locus_holdout_mask.npy"))
    sets = L.heldout_sets(idx)

    nt, nt_seeds = load_seed_ensemble(out_dir, "nt", args.require_seeds)
    ntinit, ntinit_seeds = load_seed_ensemble(out_dir, "ntinit", args.require_seeds)
    comp = np.load(os.path.join(out_dir, "LH_scores_comp.npy"))
    scores = {"nt": nt, "ntinit": ntinit, "comp": comp}

    cnn_path = args.cnn_ensemble
    if cnn_path is None:
        cnn_path = os.path.join(out_dir, "run2_frozen_20261018",
                                "LH_scores_cnn_ensemble_reconstructed.npy")
    if os.path.exists(cnn_path):
        cnn = np.load(cnn_path)
        if cnn.shape == nt.shape:
            scores["cnn"] = cnn
    print("NT seeds:", nt_seeds)
    print("random-init seeds:", ntinit_seeds)
    print("models:", sorted(scores))

    dev = pd.read_parquet("/home/vipuser/cast-v2-v0/v2_dev_labels.parquet").dropna(subset=["row"])
    piv = dev.pivot_table(index="row", columns="cell", values="label",
                          aggfunc="first").reindex(range(len(idx)))
    n_active = piv.fillna(0).sum(axis=1).values.astype(int)
    complete = piv.notna().all(axis=1).values
    chrom = idx["chrom"].astype(str).values
    mid = idx["midpoint"].values.astype(np.int64)

    metrics = []
    for cell_type, rows, labels in sets:
        for subset, subset_mask in (("unseen", holdout[rows]), ("seen", ~holdout[rows])):
            for dact, dact_mask in (
                ("all", np.ones(len(rows), dtype=bool)),
                ("dev_active", (n_active[rows] >= 1) & complete[rows]),
                ("dev_inactive", (n_active[rows] == 0) & complete[rows]),
            ):
                mask = subset_mask & dact_mask
                y = labels[mask]
                if y.sum() < 10:
                    continue
                rr = rows[mask]
                for tag, p in scores.items():
                    ap_value = average_precision_score(y, p[rr])
                    metrics.append(dict(
                        model=tag, cell=cell_type, subset=subset, dact=dact,
                        n=len(y), pos=int(y.sum()), prevalence=float(y.mean()),
                        auprc=float(ap_value), lift=float(ap_value / y.mean()),
                    ))
    metrics_df = pd.DataFrame(metrics)
    metrics_path = os.path.join(out_dir, "LH_ntinit_metrics.tsv")
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    print(metrics_df.groupby(["subset", "dact", "model"])[["auprc", "lift"]].mean().to_string())

    rng = np.random.default_rng(20261018)
    comparisons = [
        ("nt_minus_ntinit", "nt", "ntinit"),
        ("ntinit_minus_comp", "ntinit", "comp"),
        ("nt_minus_comp", "nt", "comp"),
    ]
    if "cnn" in scores:
        comparisons.append(("nt_minus_cnn", "nt", "cnn"))

    delta_rows = []
    macro_rows = []
    for subset, dact in KEY_STRATA:
        per_cell = {name: {} for name, _, _ in comparisons}
        for cell_type, rows, labels in sets:
            subset_mask = holdout[rows] if subset == "unseen" else ~holdout[rows]
            if dact == "dev_active":
                dact_mask = (n_active[rows] >= 1) & complete[rows]
            elif dact == "dev_inactive":
                dact_mask = (n_active[rows] == 0) & complete[rows]
            else:
                dact_mask = np.ones(len(rows), dtype=bool)
            mask = subset_mask & dact_mask
            y = labels[mask]
            if y.sum() < 10:
                continue
            rr = rows[mask]
            groups = block_indices(chrom[rr], mid[rr])
            for name, a_tag, b_tag in comparisons:
                point, lo, hi, skipped, draws = paired_block_delta(
                    y, scores[a_tag][rr], scores[b_tag][rr], groups, rng, args.boot
                )
                delta_rows.append(dict(
                    comparison=name, subset=subset, dact=dact, cell=cell_type,
                    n=len(y), pos=int(y.sum()), point_delta=float(point),
                    ci_lo=float(lo), ci_hi=float(hi), skipped_resamples=int(skipped),
                ))
                if np.isfinite(lo):
                    per_cell[name][cell_type] = (point, draws)

        for name, _, _ in comparisons:
            values = per_cell[name]
            if not values:
                continue
            cells = list(values)
            point = float(np.mean([values[c][0] for c in cells]))

            def macro_draw():
                picked_cells = rng.integers(0, len(cells), len(cells))
                return float(np.mean([
                    rng.choice(values[cells[i]][1]) for i in picked_cells
                ]))

            macro_samples = []
            for _ in range(args.macro_boot):
                macro_samples.append(macro_draw())
            lo = float(np.percentile(macro_samples, 2.5))
            hi = float(np.percentile(macro_samples, 97.5))
            macro_rows.append(dict(
                comparison=name, subset=subset, dact=dact, cells=len(cells),
                point_delta=point, ci_lo=lo, ci_hi=hi,
            ))

    delta_df = pd.DataFrame(delta_rows)
    macro_df = pd.DataFrame(macro_rows)
    delta_df.to_csv(os.path.join(out_dir, "LH_ntinit_deltas.tsv"), sep="\t", index=False)
    macro_df.to_csv(os.path.join(out_dir, "LH_ntinit_macro.tsv"), sep="\t", index=False)
    print("\n=== macro comparisons ===")
    print(macro_df.to_string(index=False, float_format=lambda x: "%+.4f" % x))
    print("DONE")


if __name__ == "__main__":
    main()
