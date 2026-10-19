#!/usr/bin/env python
"""Figure 1 v10.

Genome-browser hero panel, locus-by-cell activity matrix, model ladder and the
four evaluation strata. All activity values come from fig1_matrix.tsv; the
browser panel encodes FANTOM5 calls, not continuous signal.
"""
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Rectangle


HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "figures" / "main"
DATA = json.loads((HERE / "fig1_locus.json").read_text(encoding="utf-8"))
MATRIX = pd.read_csv(HERE / "fig1_matrix.tsv", sep="\t")

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 6.0,
    "axes.linewidth": 0.5,
    "savefig.dpi": 600,
})

INK = "#1B2430"
MUTED = "#6F7783"
LINE = "#D7DDE4"
NAVY = "#14508C"
STEEL = "#8FB2D4"
TEAL = "#3E8E8A"
WARM = "#C07A2C"
GRAY = "#9AA3AD"
PALE = "#EEF3F8"
BASE = {"A": "#2F855A", "C": "#2B6CB0", "G": "#B7791F", "T": "#C53030"}


def rbox(ax, x, y, w, h, fc, ec="none", lw=0.5, r=0.008, z=2, ls="-", clip=True):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={r}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z,
        linestyle=ls, clip_on=clip,
    )
    ax.add_patch(patch)
    return patch


def letter(fig, x, y, s):
    fig.text(x, y, s, fontsize=9.0, fontweight="bold", color=INK,
             ha="left", va="top")


def draw_ideogram(ax, x0, x1, yc, height=0.014):
    length = x1 - x0
    cen = x0 + 0.355 * length
    top = [(x0, yc + height), (cen - 0.014, yc + height),
           (cen, yc + 0.30 * height), (cen + 0.014, yc + height),
           (x1, yc + height)]
    bottom = [(x1, yc - height), (cen + 0.014, yc - height),
              (cen, yc - 0.30 * height), (cen - 0.014, yc - height),
              (x0, yc - height)]
    body = Polygon(top + bottom, closed=True, facecolor="#F4F6F8",
                   edgecolor="#8E979F", linewidth=0.45, zorder=2)
    ax.add_patch(body)
    bands = [(0.00, 0.06, "#DDE2E7"), (0.06, 0.15, "#F7F8FA"),
             (0.15, 0.24, "#C9D0D7"), (0.24, 0.33, "#EDF0F3"),
             (0.47, 0.56, "#DDE2E7"), (0.56, 0.66, "#C9D0D7"),
             (0.66, 0.79, "#EFF2F4"), (0.79, 0.90, "#D5DBE1"),
             (0.90, 1.00, "#C4CBD3")]
    for a, b, c in bands:
        rect = Rectangle((x0 + a * length, yc - height * 1.1),
                         (b - a) * length, height * 2.2,
                         facecolor=c, edgecolor="none", zorder=1)
        ax.add_patch(rect)
        rect.set_clip_path(body)
    body.set_zorder(3)
    body.set_facecolor("none")
    marker = x0 + 0.4074 * length
    ax.plot([marker, marker], [yc - 1.9 * height, yc + 1.9 * height],
            color="#B03A2E", lw=0.9, solid_capstyle="butt", zorder=4)
    ax.text(x0 - 0.006, yc, "chr12", fontsize=5.6, fontweight="bold",
            color=INK, ha="right", va="center")
    ax.text(marker + 0.004, yc - 2.9 * height, "54.3 Mb", fontsize=5.0,
            color="#B03A2E", ha="left", va="center")


def panel_a(fig):
    ax = fig.add_axes([0.045, 0.652, 0.925, 0.323])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    letter(fig, 0.032, 0.974, "a")

    draw_ideogram(ax, 0.690, 0.985, 0.952)

    x0, x1 = 0.158, 0.842
    span = x1 - x0
    enh_w = 606 / 2000 * span
    enh_c = 0.5 * (x0 + x1)
    enh_l, enh_r = enh_c - enh_w / 2, enh_c + enh_w / 2

    # coordinate ruler
    ax.plot([x0, x1], [0.845, 0.845], color="#B7BFC8", lw=0.7, zorder=1)
    for f in np.arange(0.05, 1.0, 0.05):
        xx = x0 + f * span
        ax.plot([xx, xx], [0.845, 0.840], color="#CDD4DB", lw=0.35, zorder=1)
    for f, lab in ((0.0, "0"), (0.25, "0.5"), (0.5, "1"), (0.75, "1.5"), (1.0, "2 kb")):
        xx = x0 + f * span
        ax.plot([xx, xx], [0.845, 0.832], color="#B7BFC8", lw=0.5)
        ax.text(xx, 0.809, lab, fontsize=5.0, color=MUTED,
                ha="center" if 0 < f < 1 else ("left" if f == 0 else "right"),
                va="top")

    # enhancer model
    rbox(ax, enh_l, 0.824, enh_w, 0.034, "#DCE8F4", NAVY, 0.5, r=0.004, z=3)
    ax.text(enh_c, 0.886, "enhancer · 606 bp", fontsize=5.2, color=NAVY,
            ha="center", va="bottom", zorder=4)
    for xx in (enh_l, enh_r):
        ax.plot([xx, xx], [0.402, 0.820], color="#F2F4F7", lw=0.5, zorder=0)

    # assay headers in the right status column
    ax.text(0.902, 0.790, "CAGE", fontsize=5.0, color=MUTED, ha="center", va="bottom")

    rows = [("K562", 0.680, 1), ("HepG2", 0.560, 0), ("GM12878", 0.440, 0)]
    for i, (cell, yy, cage) in enumerate(rows):
        ax.text(0.150, yy, cell, fontsize=5.7, color=INK, ha="right", va="center")
        ax.plot([x0, x1], [yy, yy], color="#E2E7EC", lw=0.6, zorder=1)
        rbox(ax, enh_l, yy - 0.040, enh_w, 0.080, "#BFD4E7",
             "#8FB0CE", 0.45, r=0.006, z=2)
        # CAGE call
        if cage:
            ax.add_patch(Polygon(
                [[enh_c - 0.028, yy + 0.048], [enh_c, yy + 0.112],
                 [enh_c + 0.028, yy + 0.048]],
                closed=False, fill=False, edgecolor=WARM, lw=1.0, zorder=4))
            ax.plot([enh_c - 0.030, enh_c - 0.028], [yy + 0.048, yy + 0.048],
                    color=WARM, lw=0.7)
            ax.plot([enh_c + 0.028, enh_c + 0.030], [yy + 0.048, yy + 0.048],
                    color=WARM, lw=0.7)
        else:
            ax.plot([enh_c - 0.026, enh_c + 0.026], [yy + 0.052, yy + 0.052],
                    color="#B9C0C8", lw=0.8, zorder=3)
        ax.scatter([0.902], [yy], s=13,
                   facecolor=WARM if cage else "white",
                   edgecolor=WARM if cage else "#B9C0C8",
                   linewidth=0.6, zorder=4)
        ax.text(0.919, yy, "active" if cage else "inactive", fontsize=5.0,
                color=WARM if cage else MUTED, ha="left", va="center")

    # sequence excerpt connected to the enhancer midpoint
    ax.plot([enh_c, enh_c], [0.398, 0.212], color="#C3CBD3", lw=0.5,
            linestyle=(0, (2, 2)), zorder=1)
    ax.add_patch(Circle((enh_c, 0.398), 0.0035, facecolor="white",
                        edgecolor="#C3CBD3", lw=0.5, zorder=3))
    seq = DATA["seq_enhancer"][:16]
    tile_w, gap = 0.0195, 0.0026
    total = len(seq) * tile_w + (len(seq) - 1) * gap
    sx = enh_c - total / 2
    for i, ch in enumerate(seq):
        xx = sx + i * (tile_w + gap)
        rbox(ax, xx, 0.130, tile_w, 0.052, "white", "#DCE1E7", 0.45, r=0.003, z=2)
        ax.text(xx + tile_w / 2, 0.156, ch, fontsize=5.4, color=BASE[ch],
                ha="center", va="center", fontfamily="DejaVu Sans Mono", zorder=3)
    ax.text(enh_c, 0.098, "16-bp sequence excerpt", fontsize=5.0, color=MUTED,
            ha="center", va="top")

    # compact legend for the ATAC call blocks
    rbox(ax, 0.158, 0.040, 0.030, 0.016, "#BFD4E7", "#8FB0CE", 0.45, r=0.003, z=2)
    ax.text(0.196, 0.048, "ATAC open call", fontsize=5.0, color=MUTED,
            ha="left", va="center")


ROW_GROUPS = [
    ("shared", NAVY), ("dev only", TEAL), ("held only", WARM), ("inactive", GRAY),
]


def panel_b(fig):
    cells = [c for c in MATRIX.columns
             if c not in ("enhancer_id", "set", "pattern")]
    mat = MATRIX[cells].to_numpy(int)
    ax = fig.add_axes([0.158, 0.075, 0.408, 0.445])

    n_rows, n_cols = mat.shape
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    # faint guide for the three development cells
    ax.add_patch(Rectangle((-0.5, -0.5), 3.0, n_rows,
                           facecolor="#F6F9FC", edgecolor="none", zorder=0))

    ys, xs = np.nonzero(mat)
    ax.scatter(xs, ys, s=13, color=NAVY, edgecolor="white", linewidth=0.35, zorder=3)
    iy, ix = np.nonzero(mat == 0)
    ax.scatter(ix, iy, s=3.2, color="#D6DCE3", edgecolor="none", zorder=2)

    # group rails and labels
    retained = int((MATRIX["set"] == "retained").sum())
    bounds = [0, 4, 8, 12, 16, 18, 20, 22, 24]
    for gi in range(8):
        top, bot = bounds[gi], bounds[gi + 1] - 1
        name, col = ROW_GROUPS[gi % 4]
        ax.add_patch(Rectangle((-0.72, top - 0.34), 0.10, (bot - top) + 0.68,
                               facecolor=col, edgecolor="none",
                               clip_on=False, zorder=3))
        ax.text(-0.90, (top + bot) / 2, name, fontsize=5.0, color=col,
                ha="right", va="center", clip_on=False)
    for yy in (3.5, 7.5, 11.5, 15.5, 17.5, 19.5, 21.5):
        ax.plot([-0.5, n_cols - 0.5], [yy, yy], color="#E7EBF0", lw=0.45, zorder=1)

    # locus-set braces on the right
    for top, bot, text, col in ((-0.30, retained - 0.70, "16 retained loci", NAVY),
                                (retained + 0.30, n_rows - 0.70, "8 held-out loci", WARM)):
        ax.plot([n_cols - 0.30, n_cols + 0.55, n_cols + 0.55, n_cols - 0.30],
                [top, top, bot, bot], color=col, lw=0.6, clip_on=False, zorder=4)
        ax.text(n_cols + 0.95, (top + bot) / 2, text, fontsize=5.0, color=col,
                rotation=90, rotation_mode="anchor",
                ha="center", va="center", clip_on=False)

    # cell labels
    dev = {"K562", "HepG2", "GM12878"}
    short = {"Adipocyte": "Adipo", "Endothelial": "Endo", "Bronchial": "Bronch",
             "Olfactory": "Olf", "Fibro-cardiac": "Fib-C", "Fibro-gingival": "Fib-G",
             "CD19 B": "CD19", "CD8 T": "CD8", "Skeletal muscle": "SkM",
             "Smooth muscle": "SmM", "MSC": "MSC"}
    for i, c in enumerate(cells):
        ax.text(i, -0.72, short.get(c, c), fontsize=5.0,
                color=NAVY if c in dev else TEAL,
                rotation=45, ha="left", va="bottom", clip_on=False,
                rotation_mode="anchor")

    # group headers above the rotated labels
    ax.plot([-0.45, 2.45], [-5.00, -5.00], color=NAVY, lw=0.6, clip_on=False)
    ax.text(1.0, -5.15, "development cells", fontsize=5.0, color=NAVY,
            ha="center", va="bottom", clip_on=False)
    ax.plot([2.55, n_cols - 0.45], [-5.00, -5.00], color=TEAL, lw=0.6, clip_on=False)
    ax.text((2.55 + n_cols - 0.45) / 2, -5.15, "11 held-out cell types",
            fontsize=5.0, color=TEAL, ha="center", va="bottom", clip_on=False)

    letter(fig, 0.045, 0.532, "b")
    fig.text(0.158, 0.061, "FANTOM5 loci (rows) scored in 14 primary cell types (columns)",
             fontsize=5.0, color=MUTED, ha="left", va="top")


def panel_c(fig):
    ax = fig.add_axes([0.665, 0.303, 0.315, 0.162])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    letter(fig, 0.648, 0.472, "c")

    rows = [0.795, 0.470, 0.145]

    # composition: small feature bars
    yy = rows[0]
    heights = [0.024, 0.042, 0.032, 0.050, 0.036, 0.044]
    for i, hh in enumerate(heights):
        ax.add_patch(Rectangle((0.014 + i * 0.019, yy - 0.036), 0.013, hh,
                               facecolor=STEEL if i % 2 else NAVY,
                               edgecolor="none", zorder=2))
    ax.plot([0.010, 0.128], [yy - 0.040, yy - 0.040], color=GRAY, lw=0.5)
    ax.text(0.165, yy, "Composition", fontsize=5.8, color=INK,
            fontweight="bold", ha="left", va="center")
    ax.text(0.505, yy, "20 composition features · logistic", fontsize=5.0,
            color=MUTED, ha="left", va="center")

    # CNN: stacked feature maps
    yy = rows[1]
    for i, ww in enumerate((0.116, 0.094, 0.072)):
        rbox(ax, 0.014 + (0.116 - ww) / 2, yy - 0.048 + i * 0.016, ww, 0.036,
             ["#DCE8F4", STEEL, NAVY][i], "none", r=0.004, z=2 + i)

    ax.text(0.165, yy, "DNA-CNN", fontsize=5.8, color=INK,
            fontweight="bold", ha="left", va="center")
    ax.text(0.505, yy, "3 convolution blocks", fontsize=5.0,
            color=MUTED, ha="left", va="center")

    # NT + LoRA
    yy = rows[2]
    for i in range(3):
        rbox(ax, 0.016, yy - 0.048 + i * 0.017, 0.100, 0.032,
             ["#DCE8F4", STEEL, NAVY][i], "none", r=0.004, z=2 + i)
    for i in range(2):
        rbox(ax, 0.114, yy - 0.036 + i * 0.026, 0.014, 0.018, WARM, "none",
             r=0.002, z=5)
    ax.plot([0.110, 0.110], [yy - 0.026, yy + 0.020], color=WARM, lw=0.5)

    ax.text(0.165, yy, "NT-v2 + LoRA", fontsize=5.8, color=INK,
            fontweight="bold", ha="left", va="center")
    ax.text(0.505, yy, "rank 8 · 50M backbone", fontsize=5.0,
            color=MUTED, ha="left", va="center")


def panel_d(fig):
    ax = fig.add_axes([0.665, 0.078, 0.315, 0.168])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    letter(fig, 0.648, 0.258, "d")

    ax.text(0.635, 0.985, "development-cell activity", fontsize=5.0, color=MUTED,
            ha="center", va="top")
    for xx, lab in ((0.465, "active"), (0.805, "inactive")):
        ax.text(xx, 0.870, lab, fontsize=5.2, color=INK, ha="center", va="top")

    rows = [(0.575, "unseen loci", True), (0.205, "retained loci", False)]
    for yy, lab, dashed in rows:
        ax.text(0.300, yy, lab, fontsize=5.2, color=INK, ha="right", va="center")

    for xx, active in ((0.465, True), (0.805, False)):
        for yy, _, dashed in rows:
            fc = "#E4EDF7" if active else "#FBFCFD"
            rbox(ax, xx - 0.125, yy - 0.088, 0.250, 0.176, fc,
                 "#C8D0D9", 0.5, r=0.014, z=1)
            # chromosome glyph
            bar_y = yy - 0.024
            ax.add_patch(FancyBboxPatch(
                (xx - 0.096, bar_y), 0.192, 0.048,
                boxstyle="round,pad=0,rounding_size=0.014",
                facecolor="white" if dashed else "#EEF1F4",
                edgecolor="#8E979F", linewidth=0.55,
                linestyle=(0, (1.2, 1.0)) if dashed else "-", zorder=2))
            for dx in (-0.058, 0.058):
                ax.plot([xx + dx, xx + dx], [bar_y + 0.001, bar_y + 0.047],
                        color="#C9D0D7", lw=0.4, zorder=2)
            ax.add_patch(Rectangle((xx - 0.029, bar_y + 0.004), 0.058, 0.040,
                                   facecolor=NAVY if active else "white",
                                   edgecolor=NAVY, lw=0.45, zorder=3))


def main():
    fig = plt.figure(figsize=(7.09, 4.90))
    panel_a(fig)
    panel_b(fig)
    panel_c(fig)
    panel_d(fig)

    # The manuscript embeds this artwork at 6.35 in; compensate so that the
    # smallest rendered glyph stays at or above 5 pt on the printed page.
    for artist in fig.findobj(mpl.text.Text):
        artist.set_fontsize(artist.get_fontsize() * 1.12)

    fig.savefig(OUT / "fig1_design.svg", bbox_inches="tight")
    fig.savefig(OUT / "fig1_design.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig1_design.png", bbox_inches="tight")
    fig.savefig(OUT / "fig1_design.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
