#!/usr/bin/env python
"""Sequence-neighborhood overlap audit and sensitivity analysis.

Reuses the frozen locus-holdout mask, labels and per-locus scores.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, "/home/vipuser/cast-v2-v0")
import lh_train as L  # noqa: E402

OUT = "/home/vipuser/cast-v2-locus-holdout"
BASE = "/home/vipuser/cast-v2-v0"
RNG = np.random.default_rng(20261018)


def load_scores():
    comp = np.load(os.path.join(OUT, "LH_scores_comp.npy")).astype(np.float64)
    nt = np.load(os.path.join(OUT, "LH_scores_nt_ensemble.npy")).astype(np.float64)
    cnn_files = [os.path.join(OUT, "LH_scores_cnn_seed_%d.npy" % s)
                 for s in (13, 37, 71, 101, 137)]
    cnn = np.mean([np.load(f).astype(np.float64) for f in cnn_files], axis=0)
    return {"comp": comp, "cnn": cnn, "nt": nt}


def overlap_flags(idx, hol):
    n = len(idx)
    chrom = idx["chrom"].astype(str).values
    mid = idx["midpoint"].values.astype(np.int64)
    ws = mid - 1000
    we = mid + 1000
    ov_any = np.zeros(n, dtype=bool)
    ov_max = np.zeros(n, dtype=np.int32)
    for c in np.unique(chrom):
        is_c = chrom == c
        ri = np.where(is_c & (~hol))[0]
        if len(ri) == 0:
            continue
        ri = ri[np.argsort(ws[ri])]
        rs, re = ws[ri], we[ri]
        hi = np.where(is_c & hol)[0]
        for i in hi:
            s, e = ws[i], we[i]
            k = np.where((rs < e) & (re > s))[0]
            if len(k):
                ov_any[i] = True
                ov_max[i] = int(np.maximum(0, np.minimum(re[k], e) - np.maximum(rs[k], s)).max())
    return ov_any, ov_max


def frac_for_half_recall(y, s):
    order = np.argsort(-s)
    yy = y[order]
    need = 0.5 * yy.sum()
    cum = np.cumsum(yy)
    j = int(np.searchsorted(cum, need, side="left"))
    return (j + 1) / len(y)


def metrics(y, s):
    pos = int(y.sum())
    return dict(n=len(y), pos=pos, prevalence=float(y.mean()),
                auprc=float(average_precision_score(y, s)),
                auroc=float(roc_auc_score(y, s)) if 0 < pos < len(y) else np.nan,
                frac50=float(frac_for_half_recall(y, s)))


def bootstrap_delta(rr, y, p1, p2, chrom, mid, boot=1000):
    blk = np.array([c + ":" + str(m // 1000000) for c, m in zip(chrom[rr], mid[rr])])
    ub, inv = np.unique(blk, return_inverse=True)
    groups = [np.where(inv == j)[0] for j in range(len(ub))]
    out = []
    for _ in range(boot):
        pick = RNG.integers(0, len(ub), len(ub))
        ii = np.concatenate([groups[j] for j in pick])
        if y[ii].sum() < 1:
            continue
        out.append(average_precision_score(y[ii], p1[ii]) -
                   average_precision_score(y[ii], p2[ii]))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--macro-boot", type=int, default=3000)
    ap.add_argument("--limit-cells", type=int, default=0)
    args = ap.parse_args()

    idx, _, _, _ = L.load_inputs()
    hol = np.load(os.path.join(OUT, "LH_locus_holdout_mask.npy"))
    sets = L.heldout_sets(idx)
    if args.limit_cells:
        sets = sets[:args.limit_cells]
    scores = load_scores()
    chrom = idx["chrom"].astype(str).values
    mid = idx["midpoint"].values.astype(np.int64)

    dev = pd.read_parquet(os.path.join(BASE, "v2_dev_labels.parquet")).dropna(subset=["row"])
    piv = dev.pivot_table(index="row", columns="cell", values="label", aggfunc="first").reindex(range(len(idx)))
    na = piv.fillna(0).sum(axis=1).values.astype(int)
    complete = piv.notna().all(axis=1).values

    ov_any, ov_max = overlap_flags(idx, hol)
    np.save(os.path.join(OUT, "LH_window_overlap_any.npy"), ov_any)
    np.save(os.path.join(OUT, "LH_window_overlap_max.npy"), ov_max)
    print("overlap loci: any=%d >=500=%d >=1000=%d" %
          (int(ov_any.sum()), int((ov_max >= 500).sum()), int((ov_max >= 1000).sum())), flush=True)

    # Audit table: per held-out cell, positives and overlap composition.
    audit, per_cell, delta_rows = [], [], []
    macro_cells = {}
    macro_points = {}
    for ct, rows, y in sets:
        sm = hol[rows]
        for label, keep in (
            ("base", np.ones(len(rows), dtype=bool)),
            ("disjoint_any", ~ov_any[rows]),
            ("disjoint_500", ov_max[rows] < 500),
            ("pos_only_any", ~((y == 1) & ov_any[rows])),
            ("pos_only_500", ~((y == 1) & (ov_max[rows] >= 500))),
        ):
            m = sm & keep
            yy = y[m]
            if yy.sum() < 10:
                continue
            rr = rows[m]
            rec = dict(cell=ct, definition=label, n=len(yy), pos=int(yy.sum()),
                       prevalence=float(yy.mean()))
            for tag, p in scores.items():
                mm = metrics(yy, p[rr])
                rec[tag + "_auprc"] = mm["auprc"]
                rec[tag + "_auroc"] = mm["auroc"]
                rec[tag + "_frac50"] = mm["frac50"]
            rec["nt_minus_comp"] = rec["nt_auprc"] - rec["comp_auprc"]
            per_cell.append(rec)

        # audit for unseen positives
        m = sm
        yy = y[m]
        rr = rows[m]
        p = yy == 1
        audit.append(dict(cell=ct, unseen_n=len(yy), unseen_pos=int(yy.sum()),
                          pos_ov_any=int((p & ov_any[rr]).sum()),
                          pos_ov_500=int((p & (ov_max[rr] >= 500)).sum()),
                          pos_ov_1000=int((p & (ov_max[rr] >= 1000)).sum()),
                          frac_ov_any=float((p & ov_any[rr]).sum() / max(1, p.sum())),
                          frac_ov_500=float((p & (ov_max[rr] >= 500)).sum() / max(1, p.sum()))))

    # Macro deltas and bootstrap CIs for base, disjoint_any, disjoint_500.
    for ct, rows, y in sets:
        sm = hol[rows]
        for dname, dm in (("all", np.ones(len(rows), dtype=bool)),
                          ("dev_active", (na[rows] >= 1) & complete[rows]),
                          ("dev_inactive", (na[rows] == 0) & complete[rows])):
            for label, keep in (("base", np.ones(len(rows), dtype=bool)),
                                ("disjoint_any", ~ov_any[rows]),
                                ("disjoint_500", ov_max[rows] < 500)):
                m = sm & dm & keep
                yy = y[m]
                if yy.sum() < 10:
                    continue
                rr = rows[m]
                d = bootstrap_delta(rr, yy, scores["nt"][rr], scores["comp"][rr], chrom, mid,
                                    boot=args.boot)
                if len(d) == 0:
                    continue
                key = (label, dname)
                macro_cells.setdefault(key, {})[ct] = d
                macro_points.setdefault(key, {})[ct] = float(np.mean(d))
                delta_rows.append(dict(cell=ct, definition=label, dact=dname,
                                       nt_minus_comp=float(np.mean(d)),
                                       ci_lo=float(np.percentile(d, 2.5)),
                                       ci_hi=float(np.percentile(d, 97.5)),
                                       n_boot=len(d)))
                print(ct, label, dname, "delta=%+.4f [%+.4f,%+.4f]" %
                      (np.mean(d), np.percentile(d, 2.5), np.percentile(d, 97.5)), flush=True)

    macro_rows = []
    for (label, dname), cells in macro_cells.items():
        keys = list(cells)
        samples = []
        for _ in range(args.macro_boot):
            pick = RNG.choice(len(keys), len(keys), replace=True)
            samples.append(np.mean([RNG.choice(cells[keys[i]]) for i in pick]))
        macro_rows.append(dict(definition=label, dact=dname, cells=len(keys),
                               point_delta=float(np.mean(list(macro_points[(label, dname)].values()))),
                               macro_ci_lo=float(np.percentile(samples, 2.5)),
                               macro_ci_hi=float(np.percentile(samples, 97.5))))

    pd.DataFrame(audit).to_csv(os.path.join(OUT, "LH_overlap_audit.tsv"), sep="\t", index=False)
    pd.DataFrame(per_cell).to_csv(os.path.join(OUT, "LH_overlap_metrics.tsv"), sep="\t", index=False)
    pd.DataFrame(delta_rows).to_csv(os.path.join(OUT, "LH_overlap_deltas.tsv"), sep="\t", index=False)
    pd.DataFrame(macro_rows).to_csv(os.path.join(OUT, "LH_overlap_macro.tsv"), sep="\t", index=False)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
