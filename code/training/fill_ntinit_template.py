#!/usr/bin/env python
"""Fill the random-init Results template from the final five-seed tables."""
import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/home/vipuser/cast-v2-locus-holdout")
    ap.add_argument("--template", default="ntinit_results_draft.md")
    ap.add_argument("--output", default="ntinit_results_filled.md")
    args = ap.parse_args()

    out = Path(args.out_dir)
    metrics = pd.read_csv(out / "LH_ntinit_metrics.tsv", sep="\t")
    macro = pd.read_csv(out / "LH_ntinit_macro.tsv", sep="\t")
    deltas = pd.read_csv(out / "LH_ntinit_deltas.tsv", sep="\t")

    key = (metrics.subset == "unseen") & (metrics.dact == "dev_active")
    sub = metrics[key]
    values = {}
    for model, prefix in (("comp", "COMP"), ("nt", "NT"), ("ntinit", "NTINIT")):
        values[prefix + "_AUPRC"] = sub[sub.model == model].auprc.mean()
        values[prefix + "_LIFT"] = sub[sub.model == model].lift.mean()

    def row(comparison):
        m = macro[(macro.comparison == comparison) &
                  (macro.subset == "unseen") &
                  (macro.dact == "dev_active")]
        if len(m) != 1:
            raise SystemExit("expected one macro row for %s" % comparison)
        return m.iloc[0]

    nt_minus_ntinit = row("nt_minus_ntinit")
    ntinit_minus_comp = row("ntinit_minus_comp")
    per_cell = deltas[(deltas.comparison == "nt_minus_ntinit") &
                      (deltas.subset == "unseen") &
                      (deltas.dact == "dev_active")]
    values.update({
        "NT_MINUS_NTINIT": nt_minus_ntinit.point_delta,
        "NT_MINUS_NTINIT_LO": nt_minus_ntinit.ci_lo,
        "NT_MINUS_NTINIT_HI": nt_minus_ntinit.ci_hi,
        "NTINIT_MINUS_COMP": ntinit_minus_comp.point_delta,
        "NTINIT_MINUS_COMP_LO": ntinit_minus_comp.ci_lo,
        "NTINIT_MINUS_COMP_HI": ntinit_minus_comp.ci_hi,
        "N_POSITIVE": int((per_cell.point_delta > 0).sum()),
    })

    text = Path(args.template).read_text(encoding="utf-8")
    for key_name, value in values.items():
        if isinstance(value, float):
            if "MINUS" in key_name:
                value = "%+.4f" % value
            elif "LIFT" in key_name:
                value = "%.2f" % value
            else:
                value = "%.3f" % value
        text = text.replace("<%s>" % key_name, str(value))
    Path(args.output).write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
