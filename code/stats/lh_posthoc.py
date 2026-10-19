#!/usr/bin/env python
"""Final analysis for the locus-holdout experiment.

Requires all five NT seed score files. Produces:
  LH_posthoc_metrics.tsv   per-cell metrics by subset (seen/unseen x dev activity)
  LH_posthoc_deltas.tsv    per-cell paired deltas (NT-comp, NT-CNN)
  LH_posthoc_macro.tsv     hierarchical macro delta with 95% CI
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

sys.path.insert(0, "/home/vipuser/cast-v2-v0")
import lh_train as L

OUT = "/home/vipuser/cast-v2-locus-holdout"
SEEDS = [13, 37, 71, 101, 137]


def load_scores(require_seeds=5):
    avail = [s for s in SEEDS if os.path.exists(os.path.join(OUT, "LH_scores_nt_seed_%d.npy" % s))]
    if len(avail) < require_seeds:
        print("only %d/%d NT seeds available: %s" % (len(avail), require_seeds, avail))
        sys.exit(2)
    nt = np.mean([np.load(os.path.join(OUT, "LH_scores_nt_seed_%d.npy" % s)) for s in avail], axis=0)
    out = {"nt": nt, "comp": np.load(os.path.join(OUT, "LH_scores_comp.npy"))}
    cnn_f = os.path.join(OUT, "LH_scores_cnn_ensemble.npy")
    if os.path.exists(cnn_f):
        out["cnn"] = np.load(cnn_f)
    return out, avail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--require-seeds", type=int, default=5)
    ap.add_argument("--boot", type=int, default=1000)
    args = ap.parse_args()

    idx, _, _, _ = L.load_inputs()
    hol = np.load(os.path.join(OUT, "LH_locus_holdout_mask.npy"))
    sets = L.heldout_sets(idx)
    scores, avail = load_scores(args.require_seeds)
    print("NT seeds used:", avail)

    dev = pd.read_parquet("/home/vipuser/cast-v2-v0/v2_dev_labels.parquet").dropna(subset=["row"])
    piv = dev.pivot_table(index="row", columns="cell", values="label",
                          aggfunc="first").reindex(range(len(idx)))
    na = piv.fillna(0).sum(axis=1).values.astype(int)
    complete = piv.notna().all(axis=1).values
    chrom = idx["chrom"].astype(str).values
    mid = idx["midpoint"].values.astype(np.int64)

    rows = []
    for ct, r, y in sets:
        for sname, sm in (("unseen", hol[r]), ("seen", ~hol[r])):
            for dname, dm in (("dev_active", (na[r] >= 1) & complete[r]),
                              ("dev_inactive", (na[r] == 0) & complete[r]),
                              ("all", np.ones(len(r), dtype=bool))):
                m = sm & dm
                yy = y[m]
                if yy.sum() < 10:
                    continue
                rr = r[m]
                for tag, p in scores.items():
                    a = average_precision_score(yy, p[rr])
                    rows.append(dict(model=tag, cell=ct, subset=sname, dact=dname,
                                     n=len(yy), pos=int(yy.sum()),
                                     prevalence=float(yy.mean()), auprc=float(a),
                                     lift=float(a / yy.mean())))
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(OUT, "LH_posthoc_metrics.tsv"), sep="\t", index=False)
    print(R.groupby(["subset", "dact", "model"])[["auprc", "lift"]].mean().to_string())

    # paired deltas + bootstrap on key strata
    rng = np.random.default_rng(20261018)
    deltas, macro = [], []
    for sname, dname in (("unseen", "all"), ("unseen", "dev_active"),
                         ("seen", "dev_active"), ("unseen", "dev_inactive")):
        dist = {}
        for ct, r, y in sets:
            sm = hol[r] if sname == "unseen" else ~hol[r]
            dm = (na[r] >= 1) & complete[r] if dname == "dev_active" else \
                 (na[r] == 0) & complete[r] if dname == "dev_inactive" else \
                 np.ones(len(r), dtype=bool)
            m = sm & dm
            yy = y[m]
            if yy.sum() < 10:
                continue
            rr = r[m]
            blk = np.array([c + ":" + str(p) for c, p in zip(chrom[rr], mid[rr] // 1000000)])
            ub = np.unique(blk)
            g = [np.where(blk == b)[0] for b in ub]
            d_nt_comp, d_nt_cnn = [], []
            for _ in range(args.boot):
                ii = np.concatenate([g[i] for i in rng.integers(0, len(ub), len(ub))])
                if yy[ii].sum() < 1:
                    continue
                a_nt = average_precision_score(yy[ii], scores["nt"][rr][ii])
                a_cp = average_precision_score(yy[ii], scores["comp"][rr][ii])
                d_nt_comp.append(a_nt - a_cp)
                if "cnn" in scores:
                    a_cn = average_precision_score(yy[ii], scores["cnn"][rr][ii])
                    d_nt_cnn.append(a_nt - a_cn)
            point_comp = average_precision_score(yy, scores["nt"][rr]) - average_precision_score(yy, scores["comp"][rr])
            deltas.append(dict(subset=sname, dact=dname, cell=ct,
                               nt_minus_comp=point_comp,
                               ci_lo=float(np.percentile(d_nt_comp, 2.5)),
                               ci_hi=float(np.percentile(d_nt_comp, 97.5))))
            dist[ct] = np.array(d_nt_comp)
        if dist:
            keys = list(dist)
            macro_samples = []
            for _ in range(3000):
                pick = rng.choice(len(keys), len(keys), replace=True)
                macro_samples.append(np.mean([rng.choice(dist[keys[i]]) for i in pick]))
            macro.append(dict(subset=sname, dact=dname, cells=len(keys),
                              macro_delta=float(np.mean([np.mean(dist[k]) for k in keys])),
                              ci_lo=float(np.percentile(macro_samples, 2.5)),
                              ci_hi=float(np.percentile(macro_samples, 97.5))))
    pd.DataFrame(deltas).to_csv(os.path.join(OUT, "LH_posthoc_deltas.tsv"), sep="\t", index=False)
    M = pd.DataFrame(macro)
    M.to_csv(os.path.join(OUT, "LH_posthoc_macro.tsv"), sep="\t", index=False)
    print("\n=== macro NT - comp ===")
    print(M.to_string(index=False, float_format=lambda x: "%+.4f" % x))
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
