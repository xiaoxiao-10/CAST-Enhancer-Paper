#!/usr/bin/env python
"""Recompute stronger-baseline comparisons on the label-eligible held-out sets."""
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

sys.path.insert(0, "/home/vipuser/cast-v2-v0")
import lh_train as L  # noqa: E402

MAIN = "/home/vipuser/cast-v2-locus-holdout"
BASE = "/home/vipuser/cast-v2-v0"
RNG = np.random.default_rng(20261018)


def boot(rr, yy, p1, p2, chrom, mid, B=1000):
    blk = np.array([c + ":" + str(m // 1000000) for c, m in zip(chrom, mid)])
    ub, inv = np.unique(blk, return_inverse=True)
    groups = [np.where(inv == j)[0] for j in range(len(ub))]
    out = []
    for _ in range(B):
        pick = RNG.integers(0, len(ub), len(ub))
        ii = np.concatenate([groups[j] for j in pick])
        if yy[ii].sum() < 1:
            continue
        out.append(average_precision_score(yy[ii], p1[ii]) -
                   average_precision_score(yy[ii], p2[ii]))
    return np.array(out)


def main():
    idx, _, _, _ = L.load_inputs()
    hol = np.load(os.path.join(MAIN, "LH_locus_holdout_mask.npy"))
    sets = L.heldout_sets(idx)
    chrom_all = idx["chrom"].astype(str).values
    mid_all = idx["midpoint"].values.astype(np.int64)
    scores = {
        "comp_logit": np.load(os.path.join(MAIN, "LH_scores_comp.npy")),
        "comp_xgb": np.load(os.path.join(MAIN, "LH_scores_comp_xgb.npy")),
        "motif_xgb": np.load(os.path.join(MAIN, "LH_scores_motif_xgb.npy")),
        "comb_xgb": np.load(os.path.join(MAIN, "LH_scores_comb_xgb.npy")),
        "nt": np.load(os.path.join(MAIN, "LH_scores_nt_ensemble.npy")),
    }
    dev = pd.read_parquet(os.path.join(BASE, "v2_dev_labels.parquet")).dropna(subset=["row"])
    dev["row"] = dev["row"].astype(int)
    piv = dev.pivot_table(index="row", columns="cell", values="label",
                          aggfunc="first").reindex(range(len(idx)))
    na = piv.fillna(0).sum(axis=1).values.astype(int)
    complete = piv.notna().all(axis=1).values

    rows, boot_rows = [], []
    dist = {}
    for ct, r, y in sets:
        for sname, sm in (("unseen", hol[r]), ("seen", ~hol[r])):
            for dname, dm in (("all", np.ones(len(r), dtype=bool)),
                              ("dev_active", (na[r] >= 1) & complete[r]),
                              ("dev_inactive", (na[r] == 0) & complete[r])):
                m = sm & dm
                yy = y[m]
                if yy.sum() < 10:
                    continue
                rr = r[m]
                rec = dict(cell=ct, subset=sname, dact=dname, n=len(rr), pos=int(yy.sum()))
                for name, p in scores.items():
                    rec[name] = float(average_precision_score(yy, p[rr]))
                rec["nt_minus_comp_xgb"] = rec["nt"] - rec["comp_xgb"]
                rec["nt_minus_comb_xgb"] = rec["nt"] - rec["comb_xgb"]
                rows.append(rec)
                if (sname, dname) in (("unseen", "all"), ("unseen", "dev_active"),
                                      ("seen", "dev_active"), ("unseen", "dev_inactive")):
                    d1 = boot(rr, yy, scores["nt"][rr], scores["comp_xgb"][rr],
                              chrom_all[rr], mid_all[rr])
                    d2 = boot(rr, yy, scores["nt"][rr], scores["comb_xgb"][rr],
                              chrom_all[rr], mid_all[rr])
                    if len(d1) and len(d2):
                        dist.setdefault((sname, dname), {})[ct] = (d1, d2)

    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(MAIN, "LH_correct_metrics.tsv"), sep="\t", index=False)
    M = R.groupby(["subset", "dact"]).agg(
        cells=("cell", "size"), comp_logit=("comp_logit", "mean"),
        comp_xgb=("comp_xgb", "mean"), motif_xgb=("motif_xgb", "mean"),
        comb_xgb=("comb_xgb", "mean"), nt=("nt", "mean"),
        nt_minus_comp_xgb=("nt_minus_comp_xgb", "mean"),
        nt_minus_comb_xgb=("nt_minus_comb_xgb", "mean")).reset_index()
    M.to_csv(os.path.join(MAIN, "LH_correct_macro.tsv"), sep="\t", index=False)
    print(M.to_string(index=False, float_format=lambda x: "%.4f" % x), flush=True)

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
    B.to_csv(os.path.join(MAIN, "LH_correct_boot.tsv"), sep="\t", index=False)
    print(B.to_string(index=False, float_format=lambda x: "%.4f" % x), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
