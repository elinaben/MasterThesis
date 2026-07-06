#!/usr/bin/env python3
"""Phospho-tau normalized to TOTAL tau (pTau / Tau), not to actin.

Question: does proNGF change the phosphorylated fraction of tau?

Three phospho datasets, all probed with AT8 (pSer202/pThr205) on the 680RD
channel, total tau on the 800CW channel of the same membrane:
  * Feb 2026  soma  (proNGF time course 45/60/75 min -> averaged into one proNGF bar)
  * 2026-07-02 soma
  * 2026-07-02 axon

For each lane:  pTau/Tau = (sum of phospho-tau bands) / (sum of total-tau bands).
Saves WesternBlots/figures/fig_phospho_tau_ratio.png and a values CSV.
"""
import os, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "WesternBlots", "figures"); os.makedirs(FIGDIR, exist_ok=True)
ORDER = ["Ctrl", "ProNGF", "siAPP", "siAPP+ProNGF"]

def rd(path):
    with open(os.path.join(ROOT, path), newline="") as f:
        return list(csv.DictReader(f))

def fnum(x):
    return float(x) if x not in ("", "NA", None) else 0.0

# ---- February 2026 soma -----------------------------------------------------
# phospho file = 680RD channel (pTau64 band1 + pTau56 band2); total file = 800CW
# channel (total tau band1, actin band2 -> only band1 is tau).
feb_lane = {1: "Ctrl", 2: "ProNGF", 3: "ProNGF", 4: "ProNGF", 5: "siAPP", 6: "siAPP+ProNGF"}
feb_ptau_bands, feb_tau_band1 = {}, {}
for r in rd("ImagesCSV/Feb26_p202_quantification_v2.csv"):
    feb_ptau_bands.setdefault(int(r["Lane"]), 0.0)
    feb_ptau_bands[int(r["Lane"])] += fnum(r["Area"])          # sum pTau64 + pTau56
for r in rd("ImagesCSV/Feb26_panTau_actin_quantification_v2.csv"):
    if int(r["Band"]) == 1:                                     # band1 = total tau (band2 = actin)
        feb_tau_band1[int(r["Lane"])] = fnum(r["Area"])
feb_ratio = {c: [] for c in ORDER}
for lane, cond in feb_lane.items():
    feb_ratio[cond].append(feb_ptau_bands[lane] / feb_tau_band1[lane])
feb = {c: (sum(v) / len(v)) for c, v in feb_ratio.items()}     # averages the proNGF time course

# ---- 2026-07-02 soma & axon -------------------------------------------------
def july(compartment):
    ptau, tau = {}, {}
    for r in rd("ImagesCSV/2026-07-02_ImageLab_bands.csv"):
        if r["compartment"] != compartment:
            continue
        key = r["condition"]
        if r["target"] == "phosphoTau":
            ptau[key] = ptau.get(key, 0.0) + fnum(r["adj_volume_int"])   # pTau64 + pTau56
        elif r["target"] == "totalTau":
            tau[key] = tau.get(key, 0.0) + fnum(r["adj_volume_int"])     # totalTau64 + totalTau56
    return {c: ptau.get(c, 0.0) / tau[c] for c in ORDER}

datasets = [
    ("February 2026 — soma", feb),
    ("2026-07-02 — soma", july("Soma")),
    ("2026-07-02 — axon", july("Axon")),
]

# ---- write values -----------------------------------------------------------
out_csv = os.path.join(ROOT, "WesternBlots", "phospho_tau_ratio_values.csv")
with open(out_csv, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["dataset", "condition", "pTau_over_Tau"])
    for name, d in datasets:
        for c in ORDER:
            w.writerow([name, c, round(d[c], 4)])
print("wrote", os.path.relpath(out_csv, ROOT))
for name, d in datasets:
    print(name, {c: round(d[c], 3) for c in ORDER})

# ---- figure -----------------------------------------------------------------
def color_for(cond):
    c = cond.lower().replace(" ", "")
    if "ctrl" in c:   return "#9aa0a6"
    if "siapp+" in c: return "#7a4fa3"
    if "siapp" in c:  return "#5b8c5a"
    return "#d1495b"                       # proNGF

fig, axes = plt.subplots(1, 3, figsize=(12, 4.4), sharey=False)
for ax, (title, d) in zip(axes, datasets):
    vals = [d[c] for c in ORDER]
    ax.bar(range(len(ORDER)), vals, color=[color_for(c) for c in ORDER],
           edgecolor="black", width=0.72)
    top = max(vals) if max(vals) > 0 else 1.0
    for i, v in enumerate(vals):
        ax.text(i, v + top * 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([c.replace("+", "\n+") for c in ORDER], rotation=30, ha="right", fontsize=9)
    ax.set_title(title, fontsize=11)
    ax.set_ylim(0, top * 1.20)
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("pTau / total Tau (a.u.)", fontsize=11)
fig.suptitle("Phosphorylated fraction of tau (AT8 pSer202/pThr205 ÷ total tau)",
             fontsize=13, fontweight="bold", y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(FIGDIR, "fig_phospho_tau_ratio.png")
fig.savefig(out, dpi=200, facecolor="white"); plt.close(fig)
print("wrote", os.path.relpath(out, ROOT))
