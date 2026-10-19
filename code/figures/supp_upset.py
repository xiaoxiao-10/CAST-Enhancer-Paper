#!/usr/bin/env python
"""Supplementary UpSet plot for shared and cell-specific active loci."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DATA = REPO / "data" / "derived" / "upset_positive_sets.tsv"
OUT = REPO / "figures" / "supplementary"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 6.5,
    "axes.linewidth": 0.6,
    "savefig.dpi": 600,
})

INK = "#1F2937"
MUTED = "#6B7280"
NAVY = "#0F4D92"
GRAY = "#A7ADB4"
WARM = "#B7791F"


def save(fig):
    fig.savefig(OUT / "figS_upset.svg", bbox_inches="tight")
    fig.savefig(OUT / "figS_upset.pdf", bbox_inches="tight")
    fig.savefig(OUT / "figS_upset.png", bbox_inches="tight")
    fig.savefig(OUT / "figS_upset.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main():
    df = pd.read_csv(DATA, sep="\t")
    sets = [c for c in df.columns if c != "enhancer_id"]
    mat = df[sets].to_numpy(dtype=np.int8)
    counts = {}
    for row in mat:
        bits = tuple(np.flatnonzero(row))
        if bits:
            counts[bits] = counts.get(bits, 0) + 1
    top = sorted(
        [(bits, n) for bits, n in counts.items() if len(bits) >= 2],
        key=lambda kv: kv[1], reverse=True,
    )[:25]
    set_sizes = mat.sum(axis=0)
    dev = set(range(3))
    held = set(range(3, len(sets)))

    def color(bits):
        b = set(bits)
        if b & dev and b & held:
            return NAVY
        if b & dev:
            return GRAY
        return WARM

    fig = plt.figure(figsize=(7.0, 4.7))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[1.08, 4.0], height_ratios=[1.0, 2.2],
        left=0.25, right=0.985, top=0.90, bottom=0.12,
        wspace=0.04, hspace=0.05,
    )
    ax_top = fig.add_subplot(gs[0, 1])
    ax_left = fig.add_subplot(gs[1, 0])
    ax_dot = fig.add_subplot(gs[1, 1])

    x = np.arange(len(top))
    vals = [v for _, v in top]
    ax_top.bar(x, vals, color=[color(b) for b, _ in top], width=0.72,
               edgecolor="white", linewidth=0.25)
    ax_top.set_ylabel("Intersection size", fontsize=6.5)
    ax_top.set_xlim(-0.6, len(top) - 0.4)
    ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_top.tick_params(axis="y", labelsize=5.8)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)

    ax_left.barh(np.arange(len(sets)), set_sizes, color=NAVY, height=0.62)
    ax_left.set_yticks(np.arange(len(sets)))
    ax_left.set_yticklabels(sets, fontsize=5.8)
    ax_left.invert_yaxis()
    ax_left.set_xlabel("Set size", fontsize=6.2)
    ax_left.tick_params(axis="x", labelsize=5.6)
    ax_left.spines["top"].set_visible(False)
    ax_left.spines["right"].set_visible(False)

    ax_dot.set_xlim(-0.6, len(top) - 0.4)
    ax_dot.set_ylim(len(sets) - 0.6, -0.6)
    ax_dot.set_yticks([])
    ax_dot.set_xticks([])
    ax_dot.spines["top"].set_visible(False)
    ax_dot.spines["right"].set_visible(False)
    ax_dot.spines["left"].set_visible(False)
    for xi, (bits, _) in enumerate(top):
        b = list(bits)
        if len(b) > 1:
            ax_dot.plot([xi, xi], [min(b), max(b)], color="#CBD5E1",
                        lw=0.8, zorder=1)
        ax_dot.scatter([xi] * len(b), b, s=13, color=color(bits),
                       edgecolor="white", linewidth=0.25, zorder=3)
        ax_dot.scatter([xi] * (len(sets) - len(b)),
                       [j for j in range(len(sets)) if j not in b],
                       s=4, color="#E5E7EB", zorder=0)

    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", ls="", color=NAVY,
               markersize=4, label="development + held-out"),
        Line2D([], [], marker="o", ls="", color=GRAY,
               markersize=4, label="development only"),
        Line2D([], [], marker="o", ls="", color=WARM,
               markersize=4, label="held-out only"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(0.57, 0.985), fontsize=5.9,
               frameon=False, handletextpad=0.4, columnspacing=1.0)
    fig.text(0.5, 0.025,
             "Top 25 multi-set intersections of active-locus sets; singleton sets omitted. "
             "Sets: 3 development cells and 11 held-out primary cell types.",
             ha="center", fontsize=5.8, color=MUTED)
    save(fig)


if __name__ == "__main__":
    main()
