#!/usr/bin/env python
"""Figure 2: stronger composition and motif baselines across four strata."""

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
RULE = "#D1D5DB"

MODELS = [
    ("comp_logit", "Comp.\n(logistic)", "#A7ADB4"),
    ("comp_xgb", "Comp.\n(XGBoost)", "#D97706"),
    ("motif_xgb", "Motif\n(XGBoost)", "#7C9A92"),
    ("comb_xgb", "Comp.+motif\n(XGBoost)", "#8B6BB1"),
    ("nt", "NT\n(LoRA)", "#0F4D92"),
]

PANELS = [
    ("unseen", "all", "Unseen loci, all", "a"),
    ("unseen", "dev_active", "Unseen, development-active", "b"),
    ("seen", "dev_active", "Retained, development-active", "c"),
    ("unseen", "dev_inactive", "Unseen, development-inactive", "d"),
]


def panel(ax, letter):
    ax.text(-0.16, 1.05, letter, transform=ax.transAxes,
            fontsize=8.5, fontweight="bold", color=INK,
            ha="left", va="bottom")


def style(ax):
    ax.tick_params(axis="both", labelsize=6.0, length=2, pad=2)
    ax.grid(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.7)
        ax.spines[spine].set_color(INK)


def save(fig):
    fig.savefig(OUT / "fig2_four_strata.svg", bbox_inches="tight")
    fig.savefig(OUT / "fig2_four_strata.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig2_four_strata.png", bbox_inches="tight")
    fig.savefig(OUT / "fig2_four_strata.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main():
    macro = pd.read_csv(BASE / "LH_correct_macro.tsv", sep="\t")
    boot = pd.read_csv(BASE / "LH_correct_boot.tsv", sep="\t")
    metrics = pd.read_csv(BASE / "LH_correct_metrics.tsv", sep="\t")

    fig = plt.figure(figsize=(7.25, 5.15))
    gs = fig.add_gridspec(
        2, 2, left=0.095, right=0.985, top=0.865, bottom=0.13,
        wspace=0.28, hspace=0.42,
    )
    axes = []
    for k, (subset, dact, title, letter) in enumerate(PANELS):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        axes.append(ax)
        row = macro[(macro.subset == subset) & (macro.dact == dact)].iloc[0]
        sub = metrics[(metrics.subset == subset) & (metrics.dact == dact)]
        prevalence = (sub["pos"] / sub["n"]).mean()

        vals = [row[m] for m, _, _ in MODELS]
        x = np.arange(len(MODELS))
        bars = ax.bar(x, vals, width=0.68,
                      color=[c for _, _, c in MODELS],
                      edgecolor="white", linewidth=0.45, zorder=3)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.003,
                    "%.3f" % v, ha="center", va="bottom",
                    fontsize=4.9, color=INK, rotation=0)

        ymax = max(vals) * 1.34 + prevalence
        ax.axhline(prevalence, color="#9CA3AF", lw=0.7,
                   ls=(0, (2, 2)), zorder=1)
        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.55, len(MODELS) - 0.45)
        ax.set_xticks(x)
        ax.set_xticklabels([lab for _, lab, _ in MODELS], fontsize=5.6)
        ax.set_title(title, fontsize=7.2, fontweight="bold",
                     loc="left", pad=5)
        if k % 2 == 0:
            ax.set_ylabel("Macro AUPRC", fontsize=6.8, labelpad=4)
        ax.text(0.98, 0.95,
                "n = %d cells; prevalence %.3f" % (int(row.cells), prevalence),
                transform=ax.transAxes, fontsize=5.4, color=MUTED,
                ha="right", va="top")

        b = boot[(boot.subset == subset) & (boot.dact == dact)]
        if len(b):
            b = b.iloc[0]
            ax.text(0.98, 0.86,
                    "Δ NT − comp.-XGB\n%+.4f (%+.4f, %+.4f)" % (
                        row.nt_minus_comp_xgb, b.comp_xgb_ci_lo, b.comp_xgb_ci_hi),
                    transform=ax.transAxes, fontsize=5.3, color=INK,
                    ha="right", va="top", linespacing=1.25)
        style(ax)
        panel(ax, letter)

    handles = [
        mpl.patches.Patch(facecolor=c, edgecolor="none", label=lab.replace("\n", " "))
        for _, lab, c in MODELS
    ]
    fig.legend(handles=handles, loc="upper center", ncol=5,
               bbox_to_anchor=(0.53, 0.985), fontsize=6.0,
               handlelength=1.1, handletextpad=0.35, columnspacing=1.15)
    fig.text(0.5, 0.025,
             "Bars are macro means over held-out cell types; 95% CIs for the primary contrast are shown in Fig. 3.",
             ha="center", fontsize=5.8, color=MUTED)
    save(fig)


if __name__ == "__main__":
    main()
