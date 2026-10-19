#!/usr/bin/env python
"""Preliminary locus-holdout comparison using available NT seed scores."""
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

sys.path.insert(0, "/home/vipuser/cast-v2-v0")
import lh_train as L

OUT = "/home/vipuser/cast-v2-locus-holdout"


def main():
    idx, _, _, _ = L.load_inputs()
    hol = np.load(os.path.join(OUT, "LH_locus_holdout_mask.npy"))
    sets = L.heldout_sets(idx)
    scores = {}
    for name in ("LH_scores_comp.npy", "LH_scores_cnn_ensemble.npy",
                 "LH_scores_nt_seed_13.npy", "LH_scores_nt_seed_37.npy",
                 "LH_scores_nt_ensemble.npy"):
        p = os.path.join(OUT, name)
        if os.path.exists(p) and os.path.getsize(p) > 100000:
            scores[name.replace("LH_scores_", "").replace(".npy", "")] = np.load(p)
    rows = []
    for ct, r, y in sets:
        for subset, m in (("unseen", hol[r]), ("seen", ~hol[r])):
            rr = r[m]
            yy = y[m]
            if yy.sum() < 10:
                continue
            for tag, p in scores.items():
                ap = average_precision_score(yy, p[rr])
                rows.append(dict(model=tag, cell=ct, subset=subset, n=len(yy),
                                 pos=int(yy.sum()), prevalence=float(yy.mean()),
                                 auprc=float(ap), lift=float(ap / yy.mean())))
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(OUT, "LH_prelim.tsv"), sep="\t", index=False)
    print(R.groupby(["model", "subset"])[["auprc", "lift"]].mean().to_string())
    piv = R.pivot_table(index="cell", columns=["subset", "model"], values="auprc")
    if ("unseen", "nt_seed_13") in piv.columns and ("unseen", "comp") in piv.columns:
        d = (piv[("unseen", "nt_seed_13")] - piv[("unseen", "comp")]).dropna()
        print("\nNT13 - comp, unseen loci: macro %+0.4f  (%d/%d cells positive)"
              % (d.mean(), int((d > 0).sum()), len(d)))
    if ("seen", "nt_seed_13") in piv.columns and ("seen", "comp") in piv.columns:
        d = (piv[("seen", "nt_seed_13")] - piv[("seen", "comp")]).dropna()
        print("NT13 - comp, seen loci:   macro %+0.4f  (%d/%d cells positive)"
              % (d.mean(), int((d > 0).sum()), len(d)))
    print("\nmodels found:", list(scores))


if __name__ == "__main__":
    main()
