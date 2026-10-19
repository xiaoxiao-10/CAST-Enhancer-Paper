#!/usr/bin/env python
"""NT vs stronger composition/motif baselines with paired block bootstrap."""
import os

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

MAIN = "/home/vipuser/cast-v2-locus-holdout"
BASE = "/home/vipuser/cast-v2-v0"
RNG = np.random.default_rng(20261018)


def boot(rr, y, p1, p2, chrom, mid, B=1000):
    blk = np.array([c + ":" + str(m // 1000000) for c, m in zip(chrom[rr], mid[rr])])
    ub, inv = np.unique(blk, return_inverse=True)
    groups = [np.where(inv == j)[0] for j in range(len(ub))]
    out = []
    for _ in range(B):
        pick = RNG.integers(0, len(ub), len(ub))
        ii = np.concatenate([groups[j] for j in pick])
        if y[ii].sum() < 1:
            continue
        out.append(average_precision_score(y[ii], p1[ii]) -
                   average_precision_score(y[ii], p2[ii]))
    return np.array(out)


def main():
    idx = pd.read_parquet(os.path.join(BASE, "v2_loci_index.parquet"))
    n = len(idx)
    hol = np.load(os.path.join(MAIN, "LH_locus_holdout_mask.npy"))
    chrom = idx["chrom"].astype(str).values
    mid = idx["midpoint"].values.astype(np.int64)
    scores = {
        "nt": np.load(os.path.join(MAIN, "LH_scores_nt_ensemble.npy")),
        "comp_logit": np.load(os.path.join(MAIN, "LH_scores_comp.npy")),
        "comp_xgb": np.load(os.path.join(MAIN, "LH_scores_comp_xgb.npy")),
        "motif_xgb": np.load(os.path.join(MAIN, "LH_scores_motif_xgb.npy")),
        "comb_xgb": np.load(os.path.join(MAIN, "LH_scores_comb_xgb.npy")),
    }
    dev = pd.read_parquet(os.path.join(BASE, "v2_dev_labels.parquet")).dropna(subset=["row"])
    dev["row"] = dev["row"].astype(int)
    piv = dev.pivot_table(index="row", columns="cell", values="label",
                          aggfunc="first").reindex(range(n))
    na = piv.fillna(0).sum(axis=1).values.astype(int)
    complete = piv.notna().all(axis=1).values
    labels = pd.read_csv(os.path.join(MAIN, "upset_positive_sets.tsv"), sep="\t")
    labels = labels.set_index("enhancer_id").reindex(idx.enhancer_id.values)
    held = [c for c in labels.columns if c not in ("K562", "HepG2", "GM12878")]

    rows, dist = [], {}
    for cell in held:
        y = labels[cell].fillna(0).values.astype(int)
        for sname, sm in (("unseen", hol), ("seen", ~hol)):
            for dname, dm in (("all", np.ones(n, dtype=bool)),
                              ("dev_active", (na >= 1) & complete),
                              ("dev_inactive", (na == 0) & complete)):
                m = sm & dm
                if y[m].sum() < 10:
                    continue
                rr = np.where(m)[0]
                yy = y[m]
                rec = dict(cell=cell, subset=sname, dact=dname, n=len(rr),
                           pos=int(yy.sum()))
                for name, p in scores.items():
                    rec[name] = float(average_precision_score(yy, p[m]))
                rec["nt_minus_comp_xgb"] = rec["nt"] - rec["comp_xgb"]
                rec["nt_minus_comb_xgb"] = rec["nt"] - rec["comb_xgb"]
                rows.append(rec)
                if (sname, dname) in (("unseen", "all"), ("unseen", "dev_active"),
                                      ("unseen", "dev_inactive"), ("seen", "dev_active")):
                    d1 = boot(rr, yy, scores["nt"][m], scores["comp_xgb"][m], chrom, mid)
                    d2 = boot(rr, yy, scores["nt"][m], scores["comb_xgb"][m], chrom, mid)
                    if len(d1) and len(d2):
                        dist.setdefault((sname, dname), {})[cell] = (d1, d2)

    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(MAIN, "LH_nt_vs_strong_baselines.tsv"), sep="\t", index=False)
    M = R.groupby(["subset", "dact"]).agg(
        cells=("cell", "size"),
        nt=("nt", "mean"), comp_xgb=("comp_xgb", "mean"),
        comb_xgb=("comb_xgb", "mean"),
        nt_minus_comp_xgb=("nt_minus_comp_xgb", "mean"),
        nt_minus_comb_xgb=("nt_minus_comb_xgb", "mean")).reset_index()
    M.to_csv(os.path.join(MAIN, "LH_nt_vs_strong_macro.tsv"), sep="\t", index=False)
    print(M.to_string(index=False, float_format=lambda x: "%.4f" % x), flush=True)

    boot_rows = []
    for (sname, dname), cells in dist.items():
        keys = list(cells)
        s1, s2 = [], []
        for _ in range(3000):
            pick = RNG.choice(len(keys), len(keys), replace=True)
            s1.append(np.mean([RNG.choice(cells[keys[i]][0]) for i in pick]))
            s2.append(np.mean([RNG.choice(cells[keys[i]][1]) for i in pick]))
        boot_rows.append(dict(subset=sname, dact=dname, cells=len(keys),
                              comp_xgb_ci_lo=float(np.percentile(s1, 2.5)),
                              comp_xgb_ci_hi=float(np.percentile(s1, 97.5)),
                              comb_xgb_ci_lo=float(np.percentile(s2, 2.5)),
                              comb_xgb_ci_hi=float(np.percentile(s2, 97.5))))
    B = pd.DataFrame(boot_rows)
    B.to_csv(os.path.join(MAIN, "LH_nt_vs_strong_boot.tsv"), sep="\t", index=False)
    print(B.to_string(index=False, float_format=lambda x: "%.4f" % x), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
