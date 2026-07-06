#!/usr/bin/env python3
"""Summary (averaged) phospho-tau figure: pTau / total tau, pooled across experiments.

Instead of one panel per experiment, this expresses phosphorylation as a
fold-change relative to each experiment's own control (blots are on arbitrary,
non-comparable scales), then averages the fold-changes across experiments and
plots mean +/- SEM, separately for soma and axon.

Phospho datasets (all: phospho on 680RD, total tau on 800CW, same membrane):
  * Feb 2026 #1 (2026-02-05/06)  soma + axon   -> new ImageLab reports
  * Feb 2026 #2 (2026-02-13)     soma          -> Feb26 v2 areas
  * 2026-07-02                   soma + axon   -> ImageLab bands
Note: 2026-07-02 soma control had no detectable phospho band (ratio = 0), so its
fold-change is undefined and that compartment/experiment is left out of the soma
average (flagged in the caption).
"""
import os, csv, statistics as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "WesternBlots", "figures"); os.makedirs(FIGDIR, exist_ok=True)
COND = ["Ctrl", "ProNGF", "siAPP", "siAPP+ProNGF"]

def rd(p):
    with open(os.path.join(ROOT, p), newline="") as f:
        return list(csv.DictReader(f))
def fnum(x):
    return float(x) if x not in ("", "NA", None) else 0.0
def avg(xs):
    return sum(xs) / len(xs)

# ratios[experiment][compartment][condition] = pTau/Tau
ratios = {}

# ---- Feb 2026 #1 (new reports) ----
d = {}
for r in rd("ImagesCSV/2026-02-06_phospho_totaltau_bands.csv"):
    d.setdefault(r["compartment"], {}).setdefault(r["condition"], []).append(
        float(r["phospho_adj"]) / float(r["totaltau_adj"]))
ratios["Feb2026#1"] = {comp: {c: avg(v[c]) for c in COND} for comp, v in d.items()}

# ---- Feb 2026 #2 soma (v2 area files) ----
lane_cond = {1: "Ctrl", 2: "ProNGF", 3: "ProNGF", 4: "ProNGF", 5: "siAPP", 6: "siAPP+ProNGF"}
ptau = {}
for r in rd("ImagesCSV/Feb26_p202_quantification_v2.csv"):
    ptau[int(r["Lane"])] = ptau.get(int(r["Lane"]), 0.0) + fnum(r["Area"])
tau1 = {int(r["Lane"]): fnum(r["Area"]) for r in rd("ImagesCSV/Feb26_panTau_actin_quantification_v2.csv") if int(r["Band"]) == 1}
tmp = {}
for lane, cond in lane_cond.items():
    tmp.setdefault(cond, []).append(ptau[lane] / tau1[lane])
ratios["Feb2026#2"] = {"Soma": {c: avg(tmp[c]) for c in COND}}

# ---- 2026-07-02 soma & axon ----
def july(comp):
    p, t = {}, {}
    for r in rd("ImagesCSV/2026-07-02_ImageLab_bands.csv"):
        if r["compartment"] != comp:
            continue
        if r["target"] == "phosphoTau":
            p[r["condition"]] = p.get(r["condition"], 0.0) + fnum(r["adj_volume_int"])
        elif r["target"] == "totalTau":
            t[r["condition"]] = t.get(r["condition"], 0.0) + fnum(r["adj_volume_int"])
    return {c: p.get(c, 0.0) / t[c] for c in COND}
ratios["July2026"] = {"Soma": july("Soma"), "Axon": july("Axon")}

# ---- fold-change vs each experiment's own control, then pool by compartment ----
def foldchanges(compartment):
    out = []
    for exp, comps in ratios.items():
        if compartment not in comps:
            continue
        base = comps[compartment]["Ctrl"]
        if base == 0:                       # undefined fold-change (July soma) -> skip
            print(f"  skip {exp} {compartment}: control ratio = 0 (undetectable)")
            continue
        out.append((exp, {c: comps[compartment][c] / base for c in COND}))
    return out

summary = {}
for comp in ["Soma", "Axon"]:
    fcs = foldchanges(comp)
    summary[comp] = {"n": len(fcs), "exps": [e for e, _ in fcs],
                     "mean": {c: avg([fc[c] for _, fc in fcs]) for c in COND},
                     "sem": {c: (st.pstdev([fc[c] for _, fc in fcs]) / (len(fcs) ** 0.5)) if len(fcs) > 1 else 0.0 for c in COND}}
    print(comp, "n=", summary[comp]["n"], summary[comp]["exps"])
    for c in COND:
        print(f"   {c:14s} {summary[comp]['mean'][c]:.3f} +/- {summary[comp]['sem'][c]:.3f}")

# ---- values CSV ----
with open(os.path.join(ROOT, "WesternBlots", "phospho_summary_values.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["compartment", "n_experiments", "condition", "mean_foldchange", "sem"])
    for comp in ["Soma", "Axon"]:
        for c in COND:
            w.writerow([comp, summary[comp]["n"], c, round(summary[comp]["mean"][c], 4), round(summary[comp]["sem"][c], 4)])

# ---- figure ----
x = np.arange(len(COND)); w = 0.38
colors = {"Soma": "#4c78a8", "Axon": "#e07b39"}
fig, ax = plt.subplots(figsize=(8.4, 5.0))
for i, comp in enumerate(["Soma", "Axon"]):
    m = [summary[comp]["mean"][c] for c in COND]
    e = [summary[comp]["sem"][c] for c in COND]
    ax.bar(x + (i - 0.5) * w, m, w, yerr=e, capsize=4, color=colors[comp],
           edgecolor="black", label=f"{comp} (n={summary[comp]['n']})")
ax.axhline(1.0, color="#888", lw=1, ls="--", zorder=0)
ax.set_xticks(x); ax.set_xticklabels([c.replace("+", "\n+") for c in COND], fontsize=10)
ax.set_ylabel("pTau / total tau\n(fold-change vs control)", fontsize=11)
ax.set_title("Phospho-tau summary across experiments (mean ± SEM)", fontsize=13, fontweight="bold")
ax.legend(frameon=False, fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
out = os.path.join(FIGDIR, "fig_phospho_summary.png")
fig.savefig(out, dpi=200, facecolor="white"); plt.close(fig)
print("wrote", os.path.relpath(out, ROOT))
