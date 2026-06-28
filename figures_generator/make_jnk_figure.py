#!/usr/bin/env python3
"""Generate figs/June26_JNK_figure.png — the JNK activation figure.

Panel A: p-JNK and Total JNK Western blots (stacked) with MW markers placed
close to the blots, cropped at full resolution from the original ChemiDoc scan.
Panel B: p-JNK / Total JNK ratio bar chart (Ctrl vs proNGF).

Source scan (IRDye 800CW, 2026-06-26, 1830x1465) contains two membranes:
  * p-JNK    = left  (smaller) membrane
  * Total JNK = right (larger)  membrane
In both blots: left lane = Ctrl, right lane = proNGF.
The 54/46 kDa labels mark the JNK p54/p46 isoforms (a tight doublet).

Run from anywhere:
    python3 figures_generator/make_jnk_figure.py
"""
import os
import numpy as np
from PIL import Image, ImageOps
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# repo root = parent of this script's directory
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCAN = os.path.join(ROOT, "ImagesTif", "ELINA 2026-06-26 10h55m32s(IRDye 800CW).tif")
OUT = os.path.join(ROOT, "figs", "June26_JNK_figure.png")

# full-res crop boxes (x0, y0, x1, y1) in the 1830x1465 scan
BOXES = {"pJNK": (570, 854, 864, 1052), "totalJNK": (1038, 466, 1275, 650)}

# ratio values from the June 26 quantification (ImagesCSV/June26_*_Results.csv)
GROUPS = ["Ctrl", "proNGF"]
RATIOS = [1.18, 2.60]
COLORS = ["#4878b6", "#e07b1a"]


def clean_crop(box):
    """Crop a membrane region and return it as dark-bands-on-white, contrast-stretched."""
    g = Image.open(SCAN).convert("L").crop(box)
    g = ImageOps.invert(g)                      # fluorescent bands -> dark on light
    g = ImageOps.autocontrast(g, cutoff=0.5)
    return np.asarray(g)


def band_fracs(img):
    """Return the two strongest band positions as fractions from the top of the crop."""
    dark = (255 - img).mean(axis=1).astype(float)
    dark -= dark.min()
    sm = np.convolve(dark, np.ones(5) / 5, mode="same")
    h = len(sm)
    pk = [i for i in range(4, h - 4) if sm[i] == max(sm[i - 4:i + 5]) and sm[i] > 4]
    merged = []
    for p in pk:
        if merged and p - merged[-1] < 12:
            if sm[p] > sm[merged[-1]]:
                merged[-1] = p
        else:
            merged.append(p)
    merged = sorted(merged, key=lambda p: sm[p], reverse=True)[:2]
    return sorted(p / h for p in merged)


# ---- blot images + band positions ----
p_img = clean_crop(BOXES["pJNK"])
t_img = clean_crop(BOXES["totalJNK"])
pf, tf = band_fracs(p_img), band_fracs(t_img)
P_BANDS = {"54 kDa": pf[0], "46 kDa": pf[1]}
T_BANDS = {"54 kDa": tf[0], "46 kDa": tf[1]}

# ---- figure ----
fig = plt.figure(figsize=(8.5, 9.5), dpi=200)
BL_LEFT, BL_W, BL_H = 0.30, 0.42, 0.115
P_BOT, T_BOT = 0.775, 0.600


def add_blot(img, bottom, bands, side_label, headers=False):
    ax = fig.add_axes([BL_LEFT, bottom, BL_W, BL_H])
    ax.imshow(img, cmap="gray", aspect="auto", vmin=0, vmax=255, interpolation="lanczos")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color("black")
        s.set_linewidth(1.4)

    # MW labels close to the blot; offset text so a tight doublet doesn't collide
    items = sorted(bands.items(), key=lambda kv: kv[1])      # top band first
    band_y = [bottom + BL_H * (1 - f) for _, f in items]
    min_gap = 0.026
    label_y = list(band_y)
    if len(label_y) == 2 and (label_y[0] - label_y[1]) < min_gap:
        mid = sum(label_y) / 2
        label_y = [mid + min_gap / 2, mid - min_gap / 2]
    for (name, _), by, ly in zip(items, band_y, label_y):
        fig.text(BL_LEFT - 0.020, ly, name, ha="right", va="center", fontsize=13)
        fig.add_artist(Line2D([BL_LEFT - 0.017, BL_LEFT], [ly, by], color="black", lw=1.2))

    fig.text(BL_LEFT + BL_W + 0.02, bottom + BL_H / 2, side_label,
             ha="left", va="center", fontsize=16, fontstyle="italic")
    if headers:
        fig.text(BL_LEFT + 0.28 * BL_W, bottom + BL_H + 0.012, "Ctrl",
                 ha="center", va="bottom", fontsize=14)
        fig.text(BL_LEFT + 0.72 * BL_W, bottom + BL_H + 0.012, "proNGF",
                 ha="center", va="bottom", fontsize=14)


add_blot(p_img, P_BOT, P_BANDS, "p-JNK", headers=True)
add_blot(t_img, T_BOT, T_BANDS, "Total JNK", headers=False)

fig.text(BL_LEFT + BL_W / 2, 0.945, "JNK activation: Ctrl vs proNGF",
         ha="center", va="center", fontsize=19, fontweight="bold")
fig.text(0.12, 0.945, "A", ha="center", va="center", fontsize=24, fontweight="bold")
fig.text(0.14, 0.49, "B", ha="center", va="center", fontsize=24, fontweight="bold")

axb = fig.add_axes([0.30, 0.07, 0.55, 0.40])
x = np.arange(len(GROUPS))
axb.bar(x, RATIOS, width=0.6, color=COLORS, edgecolor="black", linewidth=1.0)
for xi, r in zip(x, RATIOS):
    axb.text(xi, r + 0.05, f"{r:.2f}", ha="center", va="bottom", fontsize=14)
axb.set_xticks(x)
axb.set_xticklabels(GROUPS, fontsize=15)
axb.set_ylim(0, 3.0)
axb.set_yticks(np.arange(0, 3.01, 0.5))
axb.set_ylabel("p-JNK / Total JNK  (a.u.)", fontsize=15)
axb.yaxis.grid(True, linestyle="--", color="0.8", linewidth=1)
axb.set_axisbelow(True)
axb.spines["top"].set_visible(False)
axb.spines["right"].set_visible(False)
axb.tick_params(labelsize=12)

fig.savefig(OUT, dpi=200, facecolor="white")
print("wrote", OUT)
