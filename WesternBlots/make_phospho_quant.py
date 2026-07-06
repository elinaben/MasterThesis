#!/usr/bin/env python3
"""Quantify phospho-tau (pTau56 / actin) from two blots:
   * February 2026 (soma-only time course + siAPP)  -> Feb26 v2 Image Lab areas
   * July 2026 / 2026-07-02 (APP-knockdown, soma)    -> 2026-07-02 Image Lab volumes
Saves WesternBlots/figures/fig_phospho_quant.png
"""
import os, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "WesternBlots", "figures"); os.makedirs(FIGDIR, exist_ok=True)

def rd(path):
    with open(os.path.join(ROOT, path), newline="") as f:
        return list(csv.DictReader(f))

# ---- February 2026 (soma). Band 2 = pTau56 (p202 channel) and Band 2 = Actin (panTau+actin channel) ----
feb_lane  = {int(r["Lane"]): r["Condition"] for r in rd("ImagesCSV/LaneMapping_ExperimentJan26_2ndWB.csv")}
feb_ptau  = {int(r["Lane"]): float(r["Area"]) for r in rd("ImagesCSV/Feb26_p202_quantification_v2.csv")        if int(r["Band"]) == 2}
feb_actin = {int(r["Lane"]): float(r["Area"]) for r in rd("ImagesCSV/Feb26_panTau_actin_quantification_v2.csv") if int(r["Band"]) == 2}
feb = [(feb_lane[l], feb_ptau[l] / feb_actin[l]) for l in sorted(feb_lane)]

# ---- July 2026 / 2026-07-02 (soma) ----
jul_actin = {(r["compartment"], r["condition"]): float(r["actin"]) for r in rd("WesternBlots/2026-07-02_normalized.csv")}
jul_ptau = {}
for r in rd("ImagesCSV/2026-07-02_ImageLab_bands.csv"):
    if r["target"] == "phosphoTau" and r["band_role"] == "pTau56":
        v = r["adj_volume_int"]
        jul_ptau[(r["compartment"], r["condition"])] = float(v) if v not in ("", "NA") else 0.0
jul_order = ["Ctrl", "ProNGF", "siAPP", "siAPP+ProNGF"]
jul = [(c, jul_ptau[("Soma", c)] / jul_actin[("Soma", c)]) for c in jul_order]

def color_for(cond):
    c = cond.lower().replace(" ", "")
    if "ctrl" in c:    return "#9aa0a6"
    if "siapp+" in c:  return "#7a4fa3"
    if "siapp" in c:   return "#5b8c5a"
    return "#d1495b"                       # proNGF (any timepoint)

def panel(ax, data, title):
    labs = [c for c, _ in data]; vals = [v for _, v in data]
    ax.bar(range(len(data)), vals, color=[color_for(c) for c in labs], edgecolor="black", width=0.72)
    for i, v in enumerate(vals):
        ax.text(i, v + max(vals) * 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(data)))
    ax.set_xticklabels([c.replace("+", "\n+") for c in labs], rotation=30, ha="right", fontsize=10)
    ax.set_title(title, fontsize=12)
    ax.set_ylim(0, max(vals) * 1.18)
    ax.spines[["top", "right"]].set_visible(False)

fig, (axF, axJ) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [6, 4]})
panel(axF, feb, "February 2026 (soma)")
panel(axJ, jul, "July 2026 — 2026-07-02 (soma)")
axF.set_ylabel("pTau56 / Actin (a.u.)", fontsize=11)
fig.suptitle("Phospho-tau (pTau56) normalized to actin", fontsize=13, fontweight="bold", y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(FIGDIR, "fig_phospho_quant.png")
fig.savefig(out, dpi=200, facecolor="white"); plt.close(fig)
print("wrote", os.path.relpath(out, ROOT))
