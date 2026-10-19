#!/usr/bin/env python
"""Compare composition-only, motif-only and composition+motif XGBoost models."""
import os
import sys

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, hstack, load_npz
from sklearn.metrics import average_precision_score
from xgboost import XGBClassifier

sys.path.insert(0, "/home/vipuser/cast-v2-v0")
import lh_train as L  # noqa: E402

MAIN = "/home/vipuser/cast-v2-locus-holdout"
BASE = "/home/vipuser/cast-v2-v0"


def train_model(X, y, Xv, yv, depth, lr=0.05, n=2000):
    pos = y.sum()
    clf = XGBClassifier(
        n_estimators=n, max_depth=depth, learning_rate=lr,
        subsample=0.8, colsample_bytree=0.8, min_child_weight=1,
        scale_pos_weight=float((len(y) - pos) / max(pos, 1)),
        eval_metric="aucpr", early_stopping_rounds=50,
        tree_method="hist", n_jobs=8,
    )
    clf.fit(X, y, eval_set=[(Xv, yv)], verbose=False)
    print("depth", depth, "lr", lr, "best_iter", clf.best_iteration, flush=True)
    return clf


def main():
    idx = pd.read_parquet(os.path.join(BASE, "v2_loci_index.parquet"))
    n = len(idx)
    hol = np.load(os.path.join(MAIN, "LH_locus_holdout_mask.npy"))
    Xmotif = load_npz(os.path.join(MAIN, "motif_counts.npz")).tocsr()
    codes = np.load(os.path.join(BASE, "v2_seq_codes.npy"), mmap_mode="r")
    print("computing composition features", flush=True)
    comp = L.comp_features(codes, n)
    Xcomp = csr_matrix(comp)
    Xcomb = hstack([Xmotif, Xcomp], format="csr")

    dev = pd.read_parquet(os.path.join(BASE, "v2_dev_labels.parquet")).dropna(subset=["row"])
    dev["row"] = dev["row"].astype(int)
    dev = dev[~hol[dev.row.values]]
    tr = dev[dev.split == "train"]
    va = dev[dev.split == "validation"]
    rtr, rva = tr.row.values, va.row.values
    ytr, yva = tr.label.values, va.label.values

    models = {}
    print("training comp_xgb", flush=True)
    models["comp_xgb"] = train_model(Xcomp[rtr], ytr, Xcomp[rva], yva, depth=4)
    print("training motif_xgb", flush=True)
    models["motif_xgb"] = train_model(Xmotif[rtr], ytr, Xmotif[rva], yva, depth=6)
    print("training comb_xgb", flush=True)
    models["comb_xgb"] = train_model(Xcomb[rtr], ytr, Xcomb[rva], yva, depth=6)

    scores = {}
    for name, clf in models.items():
        X = {"comp_xgb": Xcomp, "motif_xgb": Xmotif, "comb_xgb": Xcomb}[name]
        scores[name] = clf.predict_proba(X)[:, 1]
        np.save(os.path.join(MAIN, "LH_scores_%s.npy" % name), scores[name])
    scores["comp_logit"] = np.load(os.path.join(MAIN, "LH_scores_comp.npy"))

    dev_all = pd.read_parquet(os.path.join(BASE, "v2_dev_labels.parquet")).dropna(subset=["row"])
    dev_all["row"] = dev_all["row"].astype(int)
    piv = dev_all.pivot_table(index="row", columns="cell", values="label",
                              aggfunc="first").reindex(range(n))
    na = piv.fillna(0).sum(axis=1).values.astype(int)
    complete = piv.notna().all(axis=1).values
    labels = pd.read_csv(os.path.join(MAIN, "upset_positive_sets.tsv"), sep="\t")
    labels = labels.set_index("enhancer_id").reindex(idx.enhancer_id.values)
    held = [c for c in labels.columns if c not in ("K562", "HepG2", "GM12878")]

    rows = []
    for cell in held:
        y = labels[cell].fillna(0).values.astype(int)
        for sname, sm in (("unseen", hol), ("seen", ~hol)):
            for dname, dm in (("all", np.ones(n, dtype=bool)),
                              ("dev_active", (na >= 1) & complete),
                              ("dev_inactive", (na == 0) & complete)):
                m = sm & dm
                if y[m].sum() < 10:
                    continue
                for name, s in scores.items():
                    rows.append(dict(model=name, cell=cell, subset=sname, dact=dname,
                                     n=int(m.sum()), pos=int(y[m].sum()),
                                     auprc=float(average_precision_score(y[m], s[m]))))
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(MAIN, "LH_model_compare_metrics.tsv"), sep="\t", index=False)
    M = R.groupby(["subset", "dact", "model"]).agg(
        cells=("cell", "size"), auprc=("auprc", "mean")).reset_index()
    M.to_csv(os.path.join(MAIN, "LH_model_compare_macro.tsv"), sep="\t", index=False)
    print(M.pivot_table(index=["subset", "dact"], columns="model", values="auprc").to_string(
        float_format=lambda x: "%.4f" % x), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
