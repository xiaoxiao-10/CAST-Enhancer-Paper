#!/usr/bin/env python
"""Observed (non-bootstrap) overlap-sensitivity points and secondary metrics."""
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

sys.path.insert(0, "/home/vipuser/cast-v2-v0")
sys.path.insert(0, "/home/vipuser/cast-v2-locus-holdout")
import lh_train as L  # noqa: E402
from overlap_sensitivity import load_scores, overlap_flags, metrics  # noqa: E402

OUT = "/home/vipuser/cast-v2-locus-holdout"
BASE = "/home/vipuser/cast-v2-v0"


def main():
    idx, _, _, _ = L.load_inputs()
    hol = np.load(os.path.join(OUT, "LH_locus_holdout_mask.npy"))
    sets = L.heldout_sets(idx)
    scores = load_scores()
    ov_any, ov_max = overlap_flags(idx, hol)

    dev = pd.read_parquet(os.path.join(BASE, "v2_dev_labels.parquet")).dropna(subset=["row"])
    piv = dev.pivot_table(index="row", columns="cell", values="label", aggfunc="first").reindex(range(len(idx)))
    na = piv.fillna(0).sum(axis=1).values.astype(int)
    complete = piv.notna().all(axis=1).values

    rows, secondary = [], []
    for ct, r, y in sets:
        sm = hol[r]
        for dname, dm in (("all", np.ones(len(r), dtype=bool)),
                          ("dev_active", (na[r] >= 1) & complete[r]),
                          ("dev_inactive", (na[r] == 0) & complete[r])):
            for label, keep in (("base", np.ones(len(r), dtype=bool)),
                                ("disjoint_any", ~ov_any[r]),
                                ("disjoint_500", ov_max[r] < 500)):
                m = sm & dm & keep
                yy = y[m]
                if yy.sum() < 10:
                    continue
                rr = r[m]
                rec = dict(cell=ct, definition=label, dact=dname,
                           n=len(yy), pos=int(yy.sum()))
                for tag in ("comp", "cnn", "nt"):
                    p = scores[tag][rr]
                    rec[tag + "_auprc"] = float(average_precision_score(yy, p))
                    if label == "base":
                        mm = metrics(yy, p)
                        secondary.append(dict(cell=ct, dact=dname, model=tag,
                                              n=mm["n"], pos=mm["pos"],
                                              auprc=mm["auprc"], auroc=mm["auroc"],
                                              frac50=mm["frac50"]))
                rec["nt_minus_comp"] = rec["nt_auprc"] - rec["comp_auprc"]
                rec["nt_minus_cnn"] = rec["nt_auprc"] - rec["cnn_auprc"]
                rows.append(rec)

    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(OUT, "LH_overlap_points.tsv"), sep="\t", index=False)
    M = (R.groupby(["definition", "dact"])
           .agg(cells=("cell", "size"), point_delta=("nt_minus_comp", "mean"))
           .reset_index())
    M.to_csv(os.path.join(OUT, "LH_overlap_points_macro.tsv"), sep="\t", index=False)
    pd.DataFrame(secondary).to_csv(os.path.join(OUT, "LH_secondary_metrics.tsv"), sep="\t", index=False)
    print(M.to_string(index=False, float_format=lambda x: "%+.4f" % x), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
