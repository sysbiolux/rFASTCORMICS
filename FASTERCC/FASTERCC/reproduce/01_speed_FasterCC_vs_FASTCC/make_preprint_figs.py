#!/usr/bin/env python3
"""make_preprint_figs.py -- the speed figures and tables of the preprint (2026.03.19.712885v1: Figure 3, Table 1,
Supplementary Table 1) and of the ledger HTML (Figures R1, R2, R4, S3, S5, S6, Figure 5), rebuilt from the raw data of the
full benchmark. Inputs (data/): speed_tasks.tsv (times, 1,140 tasks x 3 methods x 2 solvers), cluster_sizes.csv, and
methods_bench_tasks.csv (reactions kept per task, method and result file, from export_compare.m + compare_methods.py).
Writes into OUT (refuses to overwrite): per figure .svg/.pdf/.png + plotted data .csv + caption .md; tables .csv + .md.
'Same result' = same number of reactions kept by FasterCC and FASTCC in that task (CPLEX: results_full_cplex vs
results_full_cplex_spectra; Gurobi: same file, where the reaction IDs were also compared).
Usage: python3 make_preprint_figs.py data OUT"""
import csv, math, os, statistics as st, sys, collections
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
R = list(csv.DictReader(open(os.path.join(src, "speed_tasks.tsv")), delimiter="\t"))
SZ = {int(r["cluster"]): int(r["reactions"]) for r in csv.DictReader(open(os.path.join(src, "cluster_sizes.csv")))}
for r in R:
    for k in list(r):
        if k.endswith(("_cpx", "_grb")): r[k] = float(r[k])
    r["cluster"] = int(r["cluster"]); r["level"] = int(r["level"])
K = collections.defaultdict(dict)                          # (task name) -> {(solver, method): n}
for t in csv.DictReader(open(os.path.join(src, "methods_bench_tasks.csv"))):
    name = f"CC_c{int(t['cluster']):02d}_l{int(t['level']):03d}_r{int(t['rep']):02d}"
    if t["folder"] == "results_full_cplex" and t["method"] == "NEW": K[name][("cpx", "FC")] = int(t["n"])
    if t["folder"] == "results_full_cplex_spectra" and t["method"] in ("FASTCC", "SPECTRA"): K[name][("cpx", t["method"])] = int(t["n"])
    if t["folder"] == "results_full_gurobi_spectra": K[name][("grb", {"NEW": "FC"}.get(t["method"], t["method"]))] = int(t["n"])
CL = sorted({r["cluster"] for r in R}); LV = sorted({r["level"] for r in R})
COL = {"FasterCC": "#0072B2", "FASTCC": "#E69F00", "SPECTRA": "#009E73"}
SOL = [("cpx", "CPLEX"), ("grb", "Gurobi")]
gm = lambda v: math.exp(sum(math.log(x) for x in v) / len(v))
q = lambda v, p: st.quantiles(v, n=4)[p] if len(v) > 1 else v[0]
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e6e6", "grid.linewidth": 0.6, "axes.axisbelow": True, "svg.fonttype": "none", "pdf.fonttype": 42})
def guard(*names):
    for n in names:
        if os.path.exists(os.path.join(out, n)): raise SystemExit(n + " exists: not overwriting")
def save(fig, name, rows, header, caption):
    guard(*[f"{name}.{e}" for e in ("svg", "pdf", "png", "csv", "md")])
    for e in ("svg", "pdf", "png"): fig.savefig(os.path.join(out, f"{name}.{e}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    with open(os.path.join(out, f"{name}.csv"), "w", newline="") as fh: w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    open(os.path.join(out, f"{name}.md"), "w").write(caption.strip() + "\n"); print("written", name)
def table(name, header, rows, caption):
    guard(f"{name}.csv", f"{name}.md")
    with open(os.path.join(out, f"{name}.csv"), "w", newline="") as fh: w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    fmt = lambda x: f"{x:,}" if isinstance(x, int) and abs(x) > 999 else str(x)
    open(os.path.join(out, f"{name}.md"), "w").write(caption.strip() + "\n\n| " + " | ".join(header) + " |\n|" + "---|" * len(header) + "\n" +
                                                     "\n".join("| " + " | ".join(fmt(x) for x in r) + " |" for r in rows) + "\n"); print("written", name)
def by(cl=None, lv=None): return [r for r in R if (cl is None or r["cluster"] == cl) and (lv is None or r["level"] == lv)]
def same(task, s, other="FASTCC"):
    d = K.get(task, {}); a, b = d.get((s, "FC")), d.get((s, other))
    return None if a is None or b is None or b == 0 else a == b
# ---- Preprint Figure 3: median and IQR by number of blocked reactions, FASTCC vs FasterCC
rows = []; fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
for ax, (s, sn) in zip(axs, SOL):
    for k, mn, ls in (("FASTCC", "FASTCC", "--"), ("FC", "FasterCC", "-")):
        med = [st.median(r[f"{k}_{s}"] for r in by(lv=l)) for l in LV]; lo = [q([r[f"{k}_{s}"] for r in by(lv=l)], 0) for l in LV]
        hi = [q([r[f"{k}_{s}"] for r in by(lv=l)], 2) for l in LV]
        ax.fill_between(range(len(LV)), lo, hi, color=COL[mn], alpha=0.18, lw=0); ax.plot(range(len(LV)), med, ls, color=COL[mn], lw=1.4, marker="o", ms=3, label=f"{mn} median (IQR shaded)")
        rows += [[sn, mn, l, round(m, 2), round(a, 2), round(b, 2)] for l, m, a, b in zip(LV, med, lo, hi)]
    ax.set_xticks(range(len(LV))); ax.set_xticklabels(LV); ax.set_xlabel("blocked reactions (reactions closed at random)")
    ax.set_title(sn); ax.set_yscale("log")
axs[0].set_ylabel("time per task (s)"); axs[0].legend(frameon=False, loc="upper left")
save(fig, "figP3_runtime_by_blocked", rows, ["solver", "method", "blocked", "median_s", "q1_s", "q3_s"], """
**Figure 3 (preprint layout, full benchmark) | Runtime of FASTCC and FasterCC by number of blocked reactions.** Median time
per task (lines) and interquartile range (shaded) over the 190 tasks of each level (19 cluster networks of 21,200-183,920
reactions x 10 replicates), CPLEX (left) and Gurobi (right); log scale. Blocked reactions are created by closing 0-500
randomly chosen reactions. Data: figP3_runtime_by_blocked.csv; raw times: data/speed_tasks.tsv; code: make_preprint_figs.py.""")
# ---- Supplementary Table 1 (preprint layout): per size and level, means, SD, speed-up from means, same result
for s, sn in SOL:
    rows = []
    for c in CL:
        for l in LV:
            L = [r for r in by(c, l)]; a = [r[f"FASTCC_{s}"] for r in L]; b = [r[f"FC_{s}"] for r in L]
            sm = [same(r["task"], s) for r in L]; ok = sum(x is True for x in sm); tot = sum(x is not None for x in sm)
            rows.append([SZ[c], l, round(st.mean(a), 1), round(st.stdev(a), 1), round(st.mean(b), 1), round(st.stdev(b), 1),
                         round(st.mean(a) / st.mean(b), 1), f"{ok}/{tot}"])
    table(f"tableS1_fastcc_vs_fastercc_{sn}", ["reactions", "blocked", "FASTCC_mean_s", "FASTCC_SD_s", "FasterCC_mean_s",
          "FasterCC_SD_s", "speedup_from_means", "same_result_tasks"], rows, f"""
**Supplementary Table 1 ({sn}) | Running time of FASTCC and FasterCC by model size and number of blocked reactions.**
10 replicates per row; speed-up = mean FASTCC time / mean FasterCC time; same result = tasks in which FasterCC and FASTCC
keep the same number of reactions (Gurobi: same reaction IDs, checked in the same result file). Code: make_preprint_figs.py.""")
    if s == "cpx":
        top = sorted([r for r in rows if r[1] == 500], key=lambda r: -r[6])
        table("table1_fastcc_vs_fastercc_500blocked_CPLEX", ["reactions", "blocked", "FASTCC_mean_s", "FasterCC_mean_s", "speedup", "same_result_tasks"],
              [[r[0], r[1], r[2], r[4], r[6], r[7]] for r in top], """
**Table 1 (preprint layout, full benchmark, CPLEX) | FASTCC and FasterCC with 500 blocked reactions,** all 19 network
sizes, sorted by speed-up from the means. Code: make_preprint_figs.py.""")
# ---- Figure R1: median time per task by cluster (bars), both solvers
rows = []; fig, axs = plt.subplots(2, 1, figsize=(7.2, 4.4), sharex=True); w = 0.27
for ax, (s, sn) in zip(axs, SOL):
    for j, (k, mn) in enumerate((("FC", "FasterCC"), ("FASTCC", "FASTCC"), ("SPECTRA", "SPECTRA"))):
        y = [st.median(r[f"{k}_{s}"] for r in by(c)) for c in CL]
        ax.bar([i + (j - 1) * w for i in range(len(CL))], y, w, color=COL[mn], label=mn); rows += [[sn, mn, c, round(v, 2)] for c, v in zip(CL, y)]
    ax.set_yscale("log"); ax.set_ylabel(f"median time (s), {sn}")
axs[0].legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
axs[1].set_xticks(range(len(CL))); axs[1].set_xticklabels([f"{SZ[c]//1000}k" for c in CL], rotation=0, fontsize=7); axs[1].set_xlabel("reactions in the cluster network")
save(fig, "figR1_median_time_by_cluster", rows, ["solver", "method", "cluster", "median_s"], """
**Figure R1 | Median time per consistency-check task by cluster.** 60 tasks per cluster; CPLEX (top) and Gurobi (bottom);
log scale. Data: figR1_median_time_by_cluster.csv; code: make_preprint_figs.py.""")
# ---- Figure R2 / S5 / S6: speed-up by blocked reactions, median and IQR, vs FASTCC and vs SPECTRA
rows = []; fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
for ax, (s, sn) in zip(axs, SOL):
    for k, mn in (("FASTCC", "FASTCC"), ("SPECTRA", "SPECTRA")):
        v = [[r[f"{k}_{s}"] / r[f"FC_{s}"] for r in by(lv=l)] for l in LV]
        med = [st.median(x) for x in v]; lo = [q(x, 0) for x in v]; hi = [q(x, 2) for x in v]
        ax.fill_between(range(len(LV)), lo, hi, color=COL[mn], alpha=0.18, lw=0); ax.plot(range(len(LV)), med, "o-", ms=3, color=COL[mn], lw=1.3, label=f"{mn} / FasterCC")
        rows += [[sn, f"{mn}/FasterCC", l, round(m, 3), round(a, 3), round(b, 3)] for l, m, a, b in zip(LV, med, lo, hi)]
    ax.axhline(1, color="#666666", lw=0.8, ls="--"); ax.set_yscale("log"); ax.set_xticks(range(len(LV))); ax.set_xticklabels(LV)
    ax.set_xlabel("blocked reactions"); ax.set_title(sn)
axs[0].set_ylabel("per-task speed-up (median, IQR)"); axs[0].legend(frameon=False)
save(fig, "figR2_S5_S6_speedup_by_blocked", rows, ["solver", "ratio", "blocked", "median", "q1", "q3"], """
**Figures R2, S5, S6 | Speed-up of FasterCC by number of blocked reactions.** Per-task time of FASTCC (orange) and of
SPECTRA (green) divided by the time of FasterCC in the same task; median (line) and interquartile range (shaded) over the
190 tasks of each level; dashed line = equal time. Data: figR2_S5_S6_speedup_by_blocked.csv; code: make_preprint_figs.py.""")
# ---- Figure 5 / R4: FasterCC against SPECTRA
rows = []; fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6))
for s, sn in SOL:
    ls = "-" if s == "cpx" else "--"
    for k, mn in (("FC", "FasterCC"), ("SPECTRA", "SPECTRA")):
        y = [st.median(r[f"{k}_{s}"] for r in by(c)) for c in CL]
        axs[0].plot([SZ[c] / 1000 for c in CL], y, ls, color=COL[mn], lw=1.2, label=f"{mn} {sn}"); rows += [["A", sn, mn, c, round(v, 2)] for c, v in zip(CL, y)]
        y = [st.median(r[f"{k}_{s}"] for r in by(lv=l)) for l in LV]
        axs[2].plot(range(len(LV)), y, ls, color=COL[mn], lw=1.2, label=f"{mn} {sn}"); rows += [["C", sn, mn, l, round(v, 2)] for l, v in zip(LV, y)]
    y = [st.median(r[f"SPECTRA_{s}"] / r[f"FC_{s}"] for r in by(c)) for c in CL]
    axs[1].plot([SZ[c] / 1000 for c in CL], y, ls, color="#444444", lw=1.2, marker="o", ms=2.5, label=sn); rows += [["B", sn, "SPECTRA/FasterCC", c, round(v, 3)] for c, v in zip(CL, y)]
axs[0].set_yscale("log"); axs[0].set_xlabel("reactions (thousands)"); axs[0].set_ylabel("median time (s)"); axs[0].set_title("A  time by size", loc="left")
axs[1].axhline(1, color="#888888", lw=0.8, ls=":", label="equal time"); axs[1].set_yscale("log"); axs[1].set_xlabel("reactions (thousands)"); axs[1].set_ylabel("SPECTRA / FasterCC (median)"); axs[1].set_title("B  relative speed", loc="left"); axs[1].legend(frameon=False)
axs[2].set_yscale("log"); axs[2].set_xticks(range(len(LV))); axs[2].set_xticklabels(LV, fontsize=6); axs[2].set_xlabel("blocked reactions"); axs[2].set_title("C  time by blocked", loc="left")
axs[0].legend(frameon=False, fontsize=6); axs[1].legend(frameon=False, fontsize=6, loc="upper left"); axs[2].legend(frameon=False, fontsize=6); fig.tight_layout()
save(fig, "fig5_R4_fastercc_vs_spectra", rows, ["panel", "solver", "series", "x", "value"], """
**Figure 5 / R4 | FasterCC against SPECTRA: time and relative speed.** A, median time per task by network size; B, median
per-task ratio SPECTRA time / FasterCC time by size (above 1, FasterCC is faster); C, median time by number of blocked
reactions. Solid CPLEX, dashed Gurobi; blue FasterCC, green SPECTRA. Data: fig5_R4_fastercc_vs_spectra.csv; code: make_preprint_figs.py.""")
# ---- Figure S3: tasks with a different number of reactions kept, by cluster
rows = []; fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6), sharey=True)
for ax, (s, sn) in zip(axs, SOL):
    for j, (oth, mn) in enumerate((("FASTCC", "FASTCC"), ("SPECTRA", "SPECTRA"))):
        y = [sum(same(r["task"], s, oth) is False for r in by(c)) for c in CL]
        ax.bar([i + (j - 0.5) * 0.4 for i in range(len(CL))], y, 0.4, color=COL[mn], label=f"FasterCC vs {mn}"); rows += [[sn, mn, c, v] for c, v in zip(CL, y)]
    ax.set_xticks(range(len(CL))); ax.set_xticklabels([str(c) for c in CL], fontsize=6); ax.set_xlabel("cluster"); ax.set_title(sn)
axs[0].set_ylabel("tasks with a different count (of 60)")
from matplotlib.ticker import MaxNLocator
axs[0].yaxis.set_major_locator(MaxNLocator(integer=True))
axs[0].set_ylim(0, 3.2)
for ax in axs: ax.legend(frameon=False, fontsize=7, loc="upper left", ncol=2)
axs[1].text(0.5, 0.55, "FasterCC = FASTCC in all 1,140 tasks\nSPECTRA differs in clusters 2 and 3 only", transform=axs[1].transAxes, ha="center", fontsize=7, color="#444444")
save(fig, "figS3_tasks_with_different_result", rows, ["solver", "other_method", "cluster", "tasks_different"], """
**Figure S3 | Tasks in which the methods keep a different number of reactions,** by cluster (60 tasks each), FasterCC
against FASTCC (orange) and SPECTRA (green). With the SPEC15 reruns SPECTRA has a result in every task; every difference
in the reaction IDs is one biomass reaction (methods_bench_16Nov_diff_ids_finished.csv). Data: figS3_tasks_with_different_result.csv.""")
