#!/usr/bin/env python
"""Supplementary Fig. S3: sequence-disjoint overlap sensitivity."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "figures" / "supplementary"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 6.2,
    "axes.linewidth": 0.6,
    "savefig.dpi": 600,
})

INK = "#1B2430"
MUTED = "#6F7783"
NAVY = "#14508C"
RULE = "#D7DDE4"

ROWS = ["All held-out loci", "Remove any overlap", "Remove ≥500 bp overlap"]
DEV = [(0.0810, 0.0501, 0.1219), (0.0667, 0.0208, 0.1131), (0.0543, 0.0065, 0.0997)]
ALL = [(0.0175, 0.0101, 0.0291), (0.0184, 0.0092, 0.0324), (0.0175, 0.0084, 0.0317)]


def panel(ax, vals, title, letter, xmax, fmt="%+.3f"):
    y = np.arange(len(ROWS))[::-1]
    for yi, (p, lo, hi) in zip(y, vals):
        ax.plot([lo, hi], [yi, yi], color="#9AA7B4", lw=1.1, solid_capstyle="butt", zorder=1)
        ax.scatter([p], [yi], s=26, color=NAVY, edgecolor="white", linewidth=0.5, zorder=3)
        ax.text(hi + xmax * 0.035, yi, fmt % p, fontsize=5.8, color=NAVY,
                va="center", ha="left")
    ax.axvline(0, color="#8E979F", lw=0.6, ls=(0, (2, 2)), zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(ROWS, fontsize=6.0, color=INK)
    ax.set_ylim(-0.45, len(ROWS) - 0.55)
    ax.set_xlim(-xmax * 0.02, xmax)
    ax.set_xlabel("ΔAUPRC (NT − composition)", fontsize=6.2, labelpad=3)
    ax.tick_params(axis="x", labelsize=5.8, length=2, pad=2)
    ax.tick_params(axis="y", length=0, pad=3)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(INK)
    ax.spines["bottom"].set_color(INK)
    ax.set_title(title, fontsize=6.8, color=INK, pad=5, loc="left")
    ax.text(-0.28, 1.10, letter, transform=ax.transAxes, fontsize=8.0,
            fontweight="bold", color=INK, ha="left", va="bottom")


def main():
    fig, axes = plt.subplots(1, 2, figsize=(5.4, 2.35), sharey=True,
                             gridspec_kw={"wspace": 0.55})
    panel(axes[0], DEV, "Development-active", "a", 0.135, fmt="%+.4f")
    panel(axes[1], ALL, "All unseen loci", "b", 0.045, fmt="%+.4f")
    fig.savefig(OUT / "figS_overlap.svg", bbox_inches="tight")
    fig.savefig(OUT / "figS_overlap.pdf", bbox_inches="tight")
    fig.savefig(OUT / "figS_overlap.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
