#!/usr/bin/env python
"""Figure 3: per-cell scatter and macro forest plot against composition-XGBoost."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASE = REPO / "results" / "locus_holdout"
OUT = REPO / "figures" / "main"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.linewidth": 0.7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.0,
    "ytick.major.size": 2.0,
    "savefig.dpi": 600,
})

INK = "#1F2937"
MUTED = "#6B7280"
BLUE = "#0F4D92"
GREY = "#8B95A5"
RULE = "#D1D5DB"

ORDER = [
    ("unseen", "all", "Unseen, all", "a"),
    ("unseen", "dev_active", "Unseen, development-active", "b"),
    ("seen", "dev_active", "Retained, development-active", "c"),
    ("unseen", "dev_inactive", "Unseen, development-inactive", "d"),
]


def panel(ax, letter, x=-0.16, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=8.5, fontweight="bold", color=INK,
            ha="left", va="bottom")


def style(ax):
    ax.tick_params(axis="both", labelsize=6.2, length=2, pad=2)
    ax.grid(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.7)
        ax.spines[spine].set_color(INK)


def save(fig):
    fig.savefig(OUT / "fig3_per_cell.svg", bbox_inches="tight")
    fig.savefig(OUT / "fig3_per_cell.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig3_per_cell.png", bbox_inches="tight")
    fig.savefig(OUT / "fig3_per_cell.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main():
    metrics = pd.read_csv(BASE / "LH_correct_metrics.tsv", sep="\t")
    macro = pd.read_csv(BASE / "LH_correct_macro.tsv", sep="\t")
    boot = pd.read_csv(BASE / "LH_correct_boot.tsv", sep="\t")

    fig = plt.figure(figsize=(7.25, 3.55))
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1.08, 1.0],
        left=0.085, right=0.985, top=0.86, bottom=0.18, wspace=0.38,
    )

    # Panel a: per-cell paired scatter on the main stratum.
    ax = fig.add_subplot(gs[0, 0])
    d = metrics[(metrics.subset == "unseen") &
                (metrics.dact == "dev_active")].copy()
    x = d.comp_xgb.values
    y = d.nt.values
    pos = d["pos"].values
    better = y > x
    size = 15 + 48 * (pos - pos.min()) / max(pos.max() - pos.min(), 1)
    lo = max(0.0, min(x.min(), y.min()) * 0.9)
    hi = max(x.max(), y.max()) * 1.08
    ax.plot([lo, hi], [lo, hi], color="#94A3B8", lw=0.7,
            ls=(0, (2, 2)), zorder=1)
    ax.scatter(x[better], y[better], s=size[better], color=BLUE,
               edgecolor="white", linewidth=0.35, zorder=4,
               label="NT better")
    ax.scatter(x[~better], y[~better], s=size[~better],
               facecolors="white", edgecolors=GREY, linewidth=0.75,
               zorder=3, label="Comp.-XGB better")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Composition-XGBoost AUPRC", fontsize=6.8, labelpad=4)
    ax.set_ylabel("NT-LoRA AUPRC", fontsize=6.8, labelpad=4)
    ax.set_title("Per-cell performance on unseen\ndevelopment-active loci",
                 loc="left", fontsize=7.4, fontweight="bold", pad=6)
    ax.text(0.04, 0.95, "NT better in %d/%d cells" % (int(better.sum()), len(d)),
            transform=ax.transAxes, fontsize=6.2, color=INK,
            ha="left", va="top")
    ax.text(0.04, 0.87, "macro Δ = %+.4f" % d.nt_minus_comp_xgb.mean(),
            transform=ax.transAxes, fontsize=6.0, color=BLUE,
            ha="left", va="top")
    ax.legend(loc="lower right", fontsize=5.8, handletextpad=0.4,
              borderaxespad=0.2)
    style(ax)
    panel(ax, "a")

    # Panel b: macro forest plot across strata.
    ax = fig.add_subplot(gs[0, 1])
    ys = np.arange(len(ORDER))[::-1]
    for yv, (subset, dact, title, _) in zip(ys, ORDER):
        row = macro[(macro.subset == subset) & (macro.dact == dact)].iloc[0]
        b = boot[(boot.subset == subset) & (boot.dact == dact)]
        if len(b):
            b = b.iloc[0]
            lo_ci, hi_ci = b.comp_xgb_ci_lo, b.comp_xgb_ci_hi
        else:
            lo_ci = hi_ci = row.nt_minus_comp_xgb
        ax.errorbar(row.nt_minus_comp_xgb, yv,
                    xerr=[[row.nt_minus_comp_xgb - lo_ci],
                          [hi_ci - row.nt_minus_comp_xgb]],
                    fmt="D", color=BLUE, ecolor=INK, elinewidth=0.9,
                    capsize=2.2, markersize=4.6, zorder=4)
        ax.text(0.055, yv + 0.20, "%+.4f" % row.nt_minus_comp_xgb,
                fontsize=5.8, color=INK, ha="left", va="center")
    ax.axvline(0, color="#4B5563", lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.set_yticks(ys)
    ax.set_yticklabels([t for _, _, t, _ in ORDER], fontsize=6.0)
    ax.set_xlim(-0.065, 0.105)
    ax.set_ylim(-0.55, len(ORDER) - 0.45)
    ax.set_xlabel("ΔAUPRC (NT − composition-XGBoost)", fontsize=6.8, labelpad=5)
    ax.set_title("Macro difference by stratum", loc="left",
                 fontsize=7.4, fontweight="bold", pad=6)
    style(ax)
    panel(ax, "b", x=-0.38)

    fig.text(0.52, 0.025,
             "Diamond: macro estimate; bars: 95% paired 1 Mb block-bootstrap interval.",
             fontsize=5.6, color=MUTED, ha="center")
    save(fig)


if __name__ == "__main__":
    main()
