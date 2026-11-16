#!/usr/bin/env python3
"""make_speed_figs.py -- builds the speed figures and tables of the cluster benchmark from the raw per-task times.
Inputs (in data/):
  speed_tasks.tsv    1,140 CC tasks (clusters 2-20 x removal levels 0,10,50,100,200,500 x 10 replicates): wall time in s of
                     FasterCC (FC), FASTCC and SPECTRA, each with CPLEX (_cpx) and Gurobi (_grb); md5 f86cae10...
  cluster_sizes.csv  reactions (and reversible / irreversible) of each base cluster network, from the benchmark inputs
Writes into OUT (refuses to overwrite): for each figure a .svg, .pdf and .png (300 dpi), the data table it plots
(.csv) and its caption (.md); and table_speed_by_cluster.csv/.md (Table R1/R2 style).
Summaries used everywhere: per-task speed-up = time(other method) / time(FasterCC) in the same task and solver;
cluster value = median over its 60 tasks; overall = geometric mean of the 1,140 per-task speed-ups, and the median.
Usage: python3 make_speed_figs.py data OUT"""
import csv, math, os, statistics as st, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
R = list(csv.DictReader(open(os.path.join(src, "speed_tasks.tsv")), delimiter="\t"))
SZ = {int(r["cluster"]): int(r["reactions"]) for r in csv.DictReader(open(os.path.join(src, "cluster_sizes.csv")))}
M = [("FC", "FasterCC"), ("FASTCC", "FASTCC"), ("SPECTRA", "SPECTRA")]
COL = {"FasterCC": "#0072B2", "FASTCC": "#E69F00", "SPECTRA": "#009E73"}
SOL = [("cpx", "CPLEX"), ("grb", "Gurobi")]
for r in R:
    for k in r:
        if k.endswith(("_cpx", "_grb")): r[k] = float(r[k])
    r["cluster"] = int(r["cluster"]); r["level"] = int(r["level"])
CL = sorted({r["cluster"] for r in R}); LV = sorted({r["level"] for r in R})
gm = lambda v: math.exp(sum(math.log(x) for x in v) / len(v))
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "svg.fonttype": "none", "pdf.fonttype": 42})
def save(fig, name, rows, header, caption):
    for ext in ("svg", "pdf", "png", "csv", "md"):
        if os.path.exists(os.path.join(out, f"{name}.{ext}")): raise SystemExit(f"{name}.{ext} exists: not overwriting")
    for ext in ("svg", "pdf", "png"): fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    with open(os.path.join(out, f"{name}.csv"), "w", newline="") as fh: w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    open(os.path.join(out, f"{name}.md"), "w").write(caption.strip() + "\n"); print("written", name)
def by(cl=None, lv=None): return [r for r in R if (cl is None or r["cluster"] == cl) and (lv is None or r["level"] == lv)]
# ---------- Figure S1: median time per cluster, 3 methods, CPLEX and Gurobi
rows = []; fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
for ax, (s, sn) in zip(axs, SOL):
    for k, mn in M:
        y = [st.median(r[f"{k}_{s}"] for r in by(c)) for c in CL]
        ax.plot([SZ[c] / 1000 for c in CL], y, "o-", ms=3, lw=1.2, color=COL[mn], label=mn)
        rows += [[sn, mn, c, SZ[c], round(v, 3)] for c, v in zip(CL, y)]
    ax.set_yscale("log"); ax.set_xlabel("reactions in the cluster network (thousands)"); ax.set_title(sn)
axs[0].set_ylabel("median time per task (s)"); axs[0].legend(frameon=False)
g = {s: gm([r[f"FASTCC_{s}"] / r[f"FC_{s}"] for r in R]) for s, _ in SOL}; g2 = {s: gm([r[f"SPECTRA_{s}"] / r[f"FC_{s}"] for r in R]) for s, _ in SOL}
save(fig, "figS1_time_by_cluster", rows, ["solver", "method", "cluster", "reactions", "median_time_s"], f"""
**Figure S1 | Run time of FasterCC, FASTCC and SPECTRA on the cluster benchmark.** Median wall time per consistency-check
task (60 tasks per cluster: removal levels 0, 10, 50, 100, 200 and 500 reactions, 10 replicates each) for the 19 cluster
networks (21,200 to 183,920 reactions), with CPLEX (left) and Gurobi (right); log scale. Every method ran on the same input
files; with Gurobi all three ran in the same job, with CPLEX FasterCC ran in one job and FASTCC and SPECTRA in another;
the 15 SPECTRA tasks lost to a licence failure are the SPEC15 reruns. Over the 1,140 tasks the geometric mean of the per-task speed-up FASTCC/FasterCC is {g['cpx']:.2f}x
with CPLEX and {g['grb']:.2f}x with Gurobi; SPECTRA/FasterCC is {g2['cpx']:.2f}x and {g2['grb']:.2f}x.
Data: figS1_time_by_cluster.csv; raw times: data/speed_tasks.tsv; code: make_speed_figs.py.""")
# ---------- Figure 2: per-task speed-up by cluster (box plots)
rows = []; fig, axs = plt.subplots(2, 1, figsize=(7.2, 4.6), sharex=True)
for ax, (s, sn) in zip(axs, SOL):
    for j, (k, mn) in enumerate([("FASTCC", "FASTCC"), ("SPECTRA", "SPECTRA")]):
        data = [[r[f"{k}_{s}"] / r[f"FC_{s}"] for r in by(c)] for c in CL]
        pos = [i + (j - 0.5) * 0.36 for i in range(len(CL))]
        b = ax.boxplot(data, positions=pos, widths=0.32, patch_artist=True, showfliers=True, flierprops=dict(ms=1.5, mec=COL[mn]),
                       medianprops=dict(color="#222222", lw=1), whiskerprops=dict(color=COL[mn]), capprops=dict(color=COL[mn]))
        for p in b["boxes"]: p.set_facecolor(COL[mn]); p.set_alpha(0.55); p.set_edgecolor(COL[mn])
        ax.plot([], [], "s", color=COL[mn], alpha=0.55, label=f"{mn} / FasterCC")
        rows += [[sn, f"{mn}/FasterCC", c, round(st.median(d), 3), round(min(d), 3), round(max(d), 3)] for c, d in zip(CL, data)]
    ax.axhline(1, color="#666666", lw=0.8, ls="--"); ax.set_yscale("log"); ax.set_ylabel(f"speed-up, {sn}")
axs[0].legend(frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.0))
axs[1].set_xticks(range(len(CL))); axs[1].set_xticklabels([f"c{c}" for c in CL]); axs[1].set_xlabel("cluster network (size grows left to right)")
save(fig, "fig2_speedup_by_cluster", rows, ["solver", "ratio", "cluster", "median", "min", "max"], """
**Figure 2 | Per-task speed-up of FasterCC.** Time of FASTCC (orange) and SPECTRA (green) divided by the time of FasterCC
in the same task and with the same solver, for the 60 tasks of each cluster; top CPLEX, bottom Gurobi; log scale; the
dashed line marks equal time (above it FasterCC is faster). Boxes: quartiles and median; whiskers: 1.5 x interquartile
range; dots: tasks outside it. Data: fig2_speedup_by_cluster.csv; code: make_speed_figs.py.""")
# ---------- Figure 3: speed-up by removal level
rows = []; fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
for ax, (s, sn) in zip(axs, SOL):
    for k, mn in [("FASTCC", "FASTCC"), ("SPECTRA", "SPECTRA")]:
        y = [gm([r[f"{k}_{s}"] / r[f"FC_{s}"] for r in by(lv=l)]) for l in LV]
        ax.plot(range(len(LV)), y, "o-", ms=3, lw=1.2, color=COL[mn], label=f"{mn} / FasterCC")
        rows += [[sn, f"{mn}/FasterCC", l, round(v, 3)] for l, v in zip(LV, y)]
    ax.axhline(1, color="#666666", lw=0.8, ls="--"); ax.set_xticks(range(len(LV))); ax.set_xticklabels(LV)
    ax.set_xlabel("reactions closed at random (removal level)"); ax.set_title(sn); ax.set_yscale("log")
axs[0].set_ylabel("geometric mean speed-up (190 tasks per level)"); axs[0].legend(frameon=False)
save(fig, "fig3_speedup_by_level", rows, ["solver", "ratio", "level", "geomean_speedup"], """
**Figure 3 | Speed-up by the number of reactions closed.** Geometric mean over the 190 tasks of each removal level (19
clusters x 10 replicates) of the time of FASTCC and of SPECTRA divided by the time of FasterCC in the same task, with
CPLEX (left) and Gurobi (right); log scale; dashed line = equal time. Data: fig3_speedup_by_level.csv; code: make_speed_figs.py.""")
# ---------- Figure 4: CPLEX / Gurobi time ratio per method
rows = []; fig, ax = plt.subplots(figsize=(4.6, 2.8))
for k, mn in M:
    y = [gm([r[f"{k}_cpx"] / r[f"{k}_grb"] for r in by(c)]) for c in CL]
    ax.plot([SZ[c] / 1000 for c in CL], y, "o-", ms=3, lw=1.2, color=COL[mn], label=mn)
    rows += [[mn, c, SZ[c], round(v, 3)] for c, v in zip(CL, y)]
ax.axhline(1, color="#666666", lw=0.8, ls="--"); ax.set_yscale("log"); ax.legend(frameon=False)
ax.set_xlabel("reactions in the cluster network (thousands)"); ax.set_ylabel("time CPLEX / time Gurobi (geo. mean)")
save(fig, "fig4_solver_ratio", rows, ["method", "cluster", "reactions", "geomean_cplex_over_gurobi"], """
**Figure 4 | Effect of the LP solver.** For each method, geometric mean over the 60 tasks of each cluster of the time
with CPLEX divided by the time with Gurobi (same task); above the dashed line Gurobi is faster. The CPLEX and Gurobi runs
are separate jobs on the same inputs. Data: fig4_solver_ratio.csv; code: make_speed_figs.py.""")
# ---------- Table: per cluster medians and ratios
hdr = ["cluster", "reactions", "tasks"]; rows = []
for s, sn in SOL: hdr += [f"{sn}_FasterCC_s", f"{sn}_FASTCC_s", f"{sn}_SPECTRA_s", f"{sn}_FASTCC/FasterCC", f"{sn}_SPECTRA/FasterCC"]
for c in CL + ["all"]:
    L = R if c == "all" else by(c); row = [c, SZ.get(c, ""), len(L)]
    for s, _ in SOL:
        med = [st.median(r[f"{k}_{s}"] for r in L) for k, _ in M]
        row += [round(x, 1) for x in med] + [round(gm([r[f"FASTCC_{s}"] / r[f"FC_{s}"] for r in L]), 2),
                                           round(gm([r[f"SPECTRA_{s}"] / r[f"FC_{s}"] for r in L]), 2)]
    rows.append(row)
for ext in ("csv", "md"):
    if os.path.exists(os.path.join(out, f"table_speed_by_cluster.{ext}")): raise SystemExit("table exists: not overwriting")
with open(os.path.join(out, "table_speed_by_cluster.csv"), "w", newline="") as fh: w = csv.writer(fh); w.writerow(hdr); w.writerows(rows)
md = ["**Table | Run time per cluster.** Median time per task (s) and geometric mean of the per-task speed-ups; "
      "'all' = the 1,140 tasks. Code: make_speed_figs.py; raw data: data/speed_tasks.tsv.", "",
      "| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)] + ["| " + " | ".join(f"{x:,}" if isinstance(x, int) and x > 999 else str(x) for x in r) + " |" for r in rows]
open(os.path.join(out, "table_speed_by_cluster.md"), "w").write("\n".join(md) + "\n"); print("written table_speed_by_cluster")
