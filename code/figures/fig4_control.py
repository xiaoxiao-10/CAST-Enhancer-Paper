#!/usr/bin/env python
"""Figure 4: pretraining versus random initialization on retained and unseen loci."""
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

COL = {
    "comp": "#A7ADB4",
    "ntinit": "#42949E",
    "nt": "#0F4D92",
}
LBL = {
    "comp": "Composition",
    "ntinit": "Random-init NT",
    "nt": "Pretrained NT",
}
INK = "#1F2937"
MUTED = "#6B7280"
RULE = "#D1D5DB"

SHORT = {
    "mature adipocyte": "Mature adipocyte",
    "Endothelial Cells - Aortic": "Endothelial (aortic)",
    "Bronchial Epithelial Cell": "Bronchial epithelial",
    "Olfactory epithelial cells": "Olfactory epithelial",
    "Fibroblast - Cardiac": "Fibroblast (cardiac)",
    "Fibroblast - Gingival": "Fibroblast (gingival)",
    "CD19+ B Cells (pluriselect)": "CD19+ B cells",
    "CD8+ T Cells (pluriselect)": "CD8+ T cells",
    "Skeletal Muscle Cells": "Skeletal muscle",
    "Smooth Muscle Cells - Aortic": "Smooth muscle (aortic)",
    "Mesenchymal Stem Cells - bone marrow": "Mesenchymal stem (BM)",
}


def panel(ax, letter, x=-0.16, y=1.03):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8.5,
            fontweight="bold", color=INK, ha="left", va="bottom")


def save(fig, name):
    fig.savefig(OUT / (name + ".svg"), bbox_inches="tight")
    fig.savefig(OUT / (name + ".pdf"), bbox_inches="tight")
    fig.savefig(OUT / (name + ".png"), bbox_inches="tight")
    fig.savefig(OUT / (name + ".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def style(ax):
    ax.tick_params(axis="both", labelsize=6.5, length=2, pad=2)
    ax.grid(False)
    ax.set_facecolor("white")
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.7)
        ax.spines[spine].set_color(INK)


def main():
    metrics = pd.read_csv(BASE / "LH_ntinit_metrics.tsv", sep="\t")
    deltas = pd.read_csv(BASE / "LH_ntinit_deltas.tsv", sep="\t")
    key = (metrics.dact == "dev_active")

    means = {}
    for model in ("comp", "ntinit", "nt"):
        means[model] = {
            "seen": metrics[key & (metrics.subset == "seen") &
                          (metrics.model == model)].auprc.mean(),
            "unseen": metrics[key & (metrics.subset == "unseen") &
                            (metrics.model == model)].auprc.mean(),
        }

    fig = plt.figure(figsize=(7.0, 2.8))
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1.0, 1.12],
        wspace=0.42, left=0.085, right=0.985, top=0.86, bottom=0.20,
    )

    # Panel a: macro AUPRC slopes
    ax = fig.add_subplot(gs[0, 0])
    for model in ("ntinit", "nt"):
        sub = metrics[(metrics.dact == "dev_active") & (metrics.model == model)]
        for cell in sub.cell.unique():
            vals = []
            for subset in ("seen", "unseen"):
                vals.append(
                    sub[(sub.cell == cell) & (sub.subset == subset)].auprc.iloc[0]
                )
            ax.plot([0, 1], vals, color=COL[model], lw=0.5,
                    alpha=0.18, zorder=1)
    for model, ls in (("comp", (0, (2, 2))), ("ntinit", "-"), ("nt", "-")):
        y = [means[model]["seen"], means[model]["unseen"]]
        ax.plot([0, 1], y, color=COL[model], lw=1.6, ls=ls,
                marker="o", markersize=4.2, markeredgecolor="white",
                markeredgewidth=0.5, zorder=3)
        ax.text(1.045, y[1], "%.3f" % y[1], color=COL[model],
                fontsize=6.1, va="center", ha="left", clip_on=False)
    ax.set_xlim(-0.08, 1.30)
    ax.set_ylim(0.08, 0.30)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Retained", "Unseen"], fontsize=6.6)
    ax.set_yticks([0.10, 0.15, 0.20, 0.25, 0.30])
    ax.set_ylabel("Macro AUPRC", fontsize=6.8, labelpad=5)
    ax.set_title("Performance across locus sets", loc="left",
                 fontsize=7.4, fontweight="bold", pad=6)
    handles = [
        plt.Line2D([], [], color=COL[m], lw=1.4,
                   ls=(0, (2, 2)) if m == "comp" else "-", label=LBL[m])
        for m in ("comp", "ntinit", "nt")
    ]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.01, 0.02),
              fontsize=5.6, frameon=False, handlelength=1.4,
              handletextpad=0.4, borderaxespad=0)
    ax.text(0.04, 0.94, "thin lines: individual cell types",
            transform=ax.transAxes, fontsize=5.7, color=MUTED,
            ha="left", va="top")
    style(ax)
    panel(ax, "a")

    # Panel b: per-cell gain, retained vs unseen
    ax = fig.add_subplot(gs[0, 1])
    x = deltas[(deltas.comparison == "nt_minus_ntinit") &
               (deltas.subset == "seen") &
               (deltas.dact == "dev_active")].set_index("cell")
    y = deltas[(deltas.comparison == "nt_minus_ntinit") &
               (deltas.subset == "unseen") &
               (deltas.dact == "dev_active")].set_index("cell")
    common = [c for c in x.index if c in y.index]
    xv = x.loc[common, "point_delta"].values
    yv = y.loc[common, "point_delta"].values
    ax.plot([-0.065, 0.095], [-0.065, 0.095], color="#94A3B8",
            lw=0.7, ls=(0, (2, 2)), zorder=1)
    ax.axhline(0, color=RULE, lw=0.5, zorder=0)
    ax.axvline(0, color=RULE, lw=0.5, zorder=0)
    ax.scatter(xv, yv, s=22, color=COL["nt"], edgecolor="white",
               linewidth=0.45, zorder=3)
    for c, xx, yy in (("Fibroblast - Gingival", 0.0548, 0.2326),
                      ("Skeletal Muscle Cells", 0.0473, 0.2386),
                      ("Endothelial Cells - Aortic", 0.0210, 0.0682)):
        offsets = {
            "Fibroblast - Gingival": (5, -10),
            "Skeletal Muscle Cells": (5, 7),
            "Endothelial Cells - Aortic": (5, -2),
        }[c]
        ax.annotate(SHORT[c], (xx, yy), xytext=offsets,
                    textcoords="offset points", fontsize=5.6,
                    color=MUTED, ha="left", va="center")
    ax.text(0.03, 0.95, "above diagonal: %d/%d" % (int((yv > xv).sum()), len(common)),
            transform=ax.transAxes, ha="left", va="top",
            fontsize=6.1, color=MUTED)
    ax.set_xlim(-0.07, 0.105)
    ax.set_ylim(0.035, 0.255)
    ax.set_xticks([-0.05, 0.00, 0.05])
    ax.set_yticks([0.05, 0.10, 0.15, 0.20, 0.25])
    ax.set_xlabel("ΔAUPRC on retained loci", fontsize=6.8, labelpad=4)
    ax.set_ylabel("ΔAUPRC on unseen loci", fontsize=6.8, labelpad=5)
    ax.set_title("Per-cell gain (pretrained − random-init)", loc="left",
                 fontsize=7.4, fontweight="bold", pad=6)
    style(ax)
    panel(ax, "b")

    save(fig, "fig4_control")


if __name__ == "__main__":
    main()
