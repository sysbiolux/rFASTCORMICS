#!/usr/bin/env python3
"""plot_loopless.py -- figures of the loopless benchmark, from the per-model table written by ll_permodel2.py
(ll_15Nov_permodel.csv). Reads ONLY that CSV, so it runs anywhere (no cluster needed).
Writes into OUTDIR (refuses to overwrite):
  fig_ll_undec_by_type.(pdf|png)   undecided % per model, by model type and method (one dot per model, median bar)
  fig_ll_fastercc_stack.(pdf|png)  per model, FasterCC output: share decided by step B / closed islands / checker, and undecided
  fig_ll_fastercc_vs.(pdf|png)     per model: undecided reactions FasterCC vs FASTCC and vs SPECTRA (diagonal = equal)
Models with 0 loop-capable reactions are left out of all three (listed on screen).
Usage: python3 plot_loopless.py ll_15Nov_permodel.csv OUTDIR"""
import csv, os, sys, statistics, collections
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src, outd = sys.argv[1], sys.argv[2]; os.makedirs(outd, exist_ok=True)
R = list(csv.DictReader(open(src)))
for r in R:
    for k in ("capable", "stepB", "island", "checker", "undec"): r[k] = int(r[k])
zero = sorted({r["model"] for r in R if r["capable"] == 0}); R = [r for r in R if r["capable"] > 0]
print("left out (0 loop-capable):", ", ".join(zero) or "none")
METH = ["INPUT", "FasterCC", "FASTCC", "SPECTRA"]                                       # fixed order and colours (Okabe-Ito)
COL = {"INPUT": "#999999", "FasterCC": "#0072B2", "FASTCC": "#E69F00", "SPECTRA": "#009E73"}
TYPES = ["AGORA bacterium", "BiGG microbe", "fetal tissue", "human generic", "whole-body", "mitochondrial", "other"]
TYPES = [t for t in TYPES if any(r["type"] == t for r in R)] + sorted({r["type"] for r in R} - set(TYPES))
def save(fig, name):
    for ext in ("pdf", "png"):
        p = os.path.join(outd, f"{name}.{ext}")
        if os.path.exists(p): raise SystemExit(p + " exists: not overwriting")
        fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig); print("written", os.path.join(outd, name + ".pdf/.png"))
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e6e6", "grid.linewidth": 0.6, "axes.axisbelow": True})
# 1) undecided % per model, by type and method
fig, ax = plt.subplots(figsize=(7.2, 3.2)); w = 0.2
for i, t in enumerate(TYPES):
    for j, m in enumerate(METH):
        v = [100 * r["undec"] / r["capable"] for r in R if r["type"] == t and r["method"] == m]
        if not v: continue
        x = i + (j - 1.5) * w
        ax.scatter([x] * len(v), v, s=10, color=COL[m], alpha=0.75, linewidths=0, label=m if i == 0 else None, zorder=3)
        ax.plot([x - w * 0.4, x + w * 0.4], [statistics.median(v)] * 2, color="#333333", lw=1.2, zorder=4)
ax.set_xticks(range(len(TYPES))); ax.set_xticklabels([f"{t.replace(' ', chr(10), 1)}\n(n={len({r['model'] for r in R if r['type'] == t})})" for t in TYPES], fontsize=7)
ax.set_ylabel("undecided reactions (% of loop-capable)"); ax.set_ylim(-2, 102)
ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.28)); ax.set_title("Loop test: undecided reactions per model (bar = median)")
save(fig, "fig_ll_undec_by_type")
# 2) FasterCC output, per model: how each loop-capable reaction was decided
F = sorted([r for r in R if r["method"] == "FasterCC"], key=lambda r: (TYPES.index(r["type"]), r["model"]))
fig, ax = plt.subplots(figsize=(7.2, 0.16 * len(F) + 1.0)); y = range(len(F)); left = [0.0] * len(F)
for k, lab, c in (("stepB", "loop-free, step B (joint LP)", "#0072B2"), ("island", "loop-only, closed island", "#56B4E9"),
                  ("checker", "decided by the checker", "#009E73"), ("undec", "undecided", "#D55E00")):
    v = [100 * r[k] / r["capable"] for r in F]
    ax.barh(list(y), v, left=left, color=c, label=lab, height=0.8, edgecolor="white", linewidth=0.4)
    left = [a + b for a, b in zip(left, v)]
ax.set_yticks(list(y)); ax.set_yticklabels([f"{r['model']} ({r['capable']})" for r in F], fontsize=6); ax.invert_yaxis()
ax.set_xlabel("% of loop-capable reactions (number in brackets)"); ax.set_xlim(0, 100); ax.grid(axis="y", visible=False)
ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.6 / (0.16 * len(F) + 1.0)))
ax.set_title("FasterCC output: how each loop-capable reaction was decided")
save(fig, "fig_ll_fastercc_stack")
# 3) undecided reactions per model: FasterCC vs FASTCC, FasterCC vs SPECTRA
D = collections.defaultdict(dict)
for r in R: D[r["model"]][r["method"]] = r["undec"]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6)); fig.subplots_adjust(wspace=0.35)
for ax, other in zip(axs, ("FASTCC", "SPECTRA")):
    pts = [(d[other], d["FasterCC"]) for d in D.values() if other in d and "FasterCC" in d]
    hi = max(max(a, b) for a, b in pts) * 1.5
    ax.plot([0.8, hi], [0.8, hi], color="#888888", lw=0.8, ls="--", zorder=1, label="same number undecided")
    ax.scatter([a + 1 for a, b in pts], [b + 1 for a, b in pts], s=12, color=COL["FasterCC"], alpha=0.8, linewidths=0, zorder=3, label="one model (below the line: FasterCC leaves fewer)")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(f"undecided + 1, {other}"); ax.set_ylabel("undecided + 1, FasterCC")
    better = sum(b < a for a, b in pts); worse = sum(b > a for a, b in pts); eq = len(pts) - better - worse
    ax.set_title(f"FasterCC vs {other}\nfewer {better}, equal {eq}, more {worse} (models)", fontsize=8)
save(fig, "fig_ll_fastercc_vs")
