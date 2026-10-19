#!/usr/bin/env python
"""BOM-style motif-count baseline on the frozen locus-holdout split."""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pyfaidx
import MOODS.scan
import MOODS.tools
from scipy.sparse import csr_matrix, save_npz
from sklearn.metrics import average_precision_score
from xgboost import XGBClassifier

MAIN = "/home/vipuser/cast-v2-locus-holdout"
BASE = "/home/vipuser/cast-v2-v0"
MEME = os.path.join(MAIN, "jaspar2024_vertebrates.meme")
FASTA = os.path.join(BASE, "hg38.fa")
RNG = np.random.default_rng(20261018)


def read_meme(path):
    mats, names, cur = [], [], None
    for line in open(path):
        if line.startswith("MOTIF"):
            if cur is not None and cur:
                mats.append([[cur[i][b] for i in range(len(cur))] for b in range(4)])
                names.append(name)
            parts = line.split()
            name = parts[2] if len(parts) > 2 else parts[1]
            cur = []
        elif cur is not None and line.strip():
            parts = line.split()
            if len(parts) == 4:
                try:
                    cur.append([float(x) for x in parts])
                except ValueError:
                    pass
    if cur:
        mats.append([[cur[i][b] for i in range(len(cur))] for b in range(4)])
        names.append(name)
    return mats, names


def scan_counts(args):
    idx = pd.read_parquet(os.path.join(BASE, "v2_loci_index.parquet"))
    n = len(idx) if not args.limit else min(args.limit, len(idx))
    mats, names = read_meme(MEME)
    if args.motifs:
        mats, names = mats[:args.motifs], names[:args.motifs]
    bg = MOODS.tools.flat_bg(4)
    th = [MOODS.tools.threshold_from_p(x, bg, args.pvalue) for x in mats]
    sc = MOODS.scan.Scanner(max(len(x[0]) for x in mats))
    sc.set_motifs(mats, bg, th)
    fa = pyfaidx.Fasta(FASTA)

    chrom = idx["chrom"].astype(str).values
    mid = idx["midpoint"].values.astype(np.int64)
    rows, cols, vals = [], [], []
    for i in range(n):
        seq = str(fa[chrom[i]][mid[i] - 1000:mid[i] + 1000]).upper().replace("N", "A")
        hits = sc.scan(seq)
        for j, h in enumerate(hits):
            if h:
                rows.append(i)
                cols.append(j)
                vals.append(len(h))
        if i and i % 2000 == 0:
            print("scanned", i, flush=True)
    X = csr_matrix((vals, (rows, cols)), shape=(n, len(mats)))
    save_npz(os.path.join(MAIN, "motif_counts.npz"), X)
    np.save(os.path.join(MAIN, "motif_count_cols.npy"), np.array([len(vals), X.nnz]))
    with open(os.path.join(MAIN, "motif_names.txt"), "w") as fh:
        fh.write("\n".join(names))
    print("matrix", X.shape, "nnz", X.nnz, flush=True)
    return idx.iloc[:n].copy(), X


def evaluate(args, idx, X):
    hol = np.load(os.path.join(MAIN, "LH_locus_holdout_mask.npy"))[:len(idx)]
    dev = pd.read_parquet(os.path.join(BASE, "v2_dev_labels.parquet")).dropna(subset=["row"])
    dev["row"] = dev["row"].astype(int)
    dev = dev[~hol[dev.row.values]]
    tr = dev[dev.split == "train"]
    va = dev[dev.split == "validation"]
    pos = tr.label.sum()
    clf = XGBClassifier(
        n_estimators=args.n_estimators, max_depth=4, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        scale_pos_weight=float((len(tr) - pos) / max(pos, 1)),
        eval_metric="aucpr", early_stopping_rounds=30,
        tree_method="hist", n_jobs=8,
    )
    clf.fit(X[tr.row.values], tr.label.values,
            eval_set=[(X[va.row.values], va.label.values)], verbose=False)
    scores = clf.predict_proba(X)[:, 1]
    np.save(os.path.join(MAIN, "LH_scores_motif.npy"), scores)
    print("model done, best_iter", clf.best_iteration, flush=True)

    piv = dev.pivot_table(index="row", columns="cell", values="label",
                          aggfunc="first").reindex(range(len(idx)))
    # use all development labels (including validation/test) for activity strata
    dev_all = pd.read_parquet(os.path.join(BASE, "v2_dev_labels.parquet")).dropna(subset=["row"])
    dev_all["row"] = dev_all["row"].astype(int)
    piv = dev_all.pivot_table(index="row", columns="cell", values="label",
                              aggfunc="first").reindex(range(len(idx)))
    na = piv.fillna(0).sum(axis=1).values.astype(int)
    complete = piv.notna().all(axis=1).values

    labels = pd.read_csv(os.path.join(MAIN, "upset_positive_sets.tsv"), sep="\t")
    labels = labels.set_index("enhancer_id").reindex(idx.enhancer_id.values)
    held = [c for c in labels.columns if c not in ("K562", "HepG2", "GM12878")]
    comp = np.load(os.path.join(MAIN, "LH_scores_comp.npy"))[:len(idx)]

    rows = []
    for cell in held:
        y = labels[cell].fillna(0).values.astype(int)
        for sname, sm in (("unseen", hol), ("seen", ~hol)):
            for dname, dm in (("all", np.ones(len(idx), dtype=bool)),
                              ("dev_active", (na >= 1) & complete),
                              ("dev_inactive", (na == 0) & complete)):
                m = sm & dm
                if y[m].sum() < 10:
                    continue
                rows.append(dict(
                    cell=cell, subset=sname, dact=dname, n=int(m.sum()),
                    pos=int(y[m].sum()),
                    comp_auprc=float(average_precision_score(y[m], comp[m])),
                    motif_auprc=float(average_precision_score(y[m], scores[m])),
                ))
    R = pd.DataFrame(rows)
    R["delta"] = R.motif_auprc - R.comp_auprc
    R.to_csv(os.path.join(MAIN, "LH_motif_metrics.tsv"), sep="\t", index=False)
    M = R.groupby(["subset", "dact"]).agg(
        cells=("cell", "size"), comp_auprc=("comp_auprc", "mean"),
        motif_auprc=("motif_auprc", "mean"), delta=("delta", "mean")).reset_index()
    M.to_csv(os.path.join(MAIN, "LH_motif_macro.tsv"), sep="\t", index=False)
    print(M.to_string(index=False, float_format=lambda x: "%.4f" % x), flush=True)
    print("DONE", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--motifs", type=int, default=0)
    ap.add_argument("--pvalue", type=float, default=1e-4)
    ap.add_argument("--n-estimators", type=int, default=600)
    args = ap.parse_args()
    idx, X = scan_counts(args)
    if not args.limit:
        evaluate(args, idx, X)


if __name__ == "__main__":
    main()
