#!/usr/bin/env python3
"""make_spectra_figs.py -- FasterCC compared with SPECTRA (and FASTCC): speed per task, speed by network size and number
of blocked reactions, real-model run times, agreement with FVA, and reversible reactions left in the output.
Inputs (data/):
  speed_tasks.tsv          time per task (1,140 cluster tasks x 3 methods x 2 solvers; SPECTRA times of the 15 licence-
                           failed tasks are those of the SPEC15 reruns)
  cluster_sizes.csv        reactions per cluster network
  methods_bench_tasks.csv  reactions kept, reversible counts, per task / method / result file (export_compare.m +
                           compare_methods.py)
  spec15_tasks.tsv         the 15 SPEC15 reruns (export_compare.m 'bench' on results_spec15_*; copied from the cluster
                           output of 16 Nov 03:46)
  methods_real_16Nov_models.csv  58 real models, the three methods (results_real/{NEW,FASTCC,SPECTRA})
Writes into OUT (refuses to overwrite): per figure .svg/.pdf/.png + plotted data .csv + caption .md; tables .csv + .md.
Usage: python3 make_spectra_figs.py data OUT"""
import csv, math, os, statistics as st, sys, collections
import matplotlib, matplotlib.ticker; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
rd = lambda f, d=",": list(csv.DictReader(open(os.path.join(src, f)), delimiter=d))
R = rd("speed_tasks.tsv", "\t"); SZ = {int(r["cluster"]): int(r["reactions"]) for r in rd("cluster_sizes.csv")}
for r in R:
    for k in list(r):
        if k.endswith(("_cpx", "_grb")): r[k] = float(r[k])
    r["cluster"] = int(r["cluster"]); r["level"] = int(r["level"])
BT = rd("methods_bench_tasks.csv"); S15 = rd("spec15_tasks.tsv", "\t"); RM = rd("methods_real_16Nov_models.csv")
CL = sorted({r["cluster"] for r in R}); LV = sorted({r["level"] for r in R})
COL = {"FasterCC": "#0072B2", "FASTCC": "#E69F00", "SPECTRA": "#009E73"}; MK = {"FasterCC": "o", "FASTCC": "s", "SPECTRA": "^"}
SOL = [("cpx", "CPLEX"), ("grb", "Gurobi")]
gm = lambda v: math.exp(sum(math.log(x) for x in v) / len(v))
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

# ---------------- SP1: per-task time, FasterCC vs SPECTRA (both solvers)
rows = []; fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.4), sharex=True, sharey=True); norm = LogNorm(min(SZ.values()), max(SZ.values()))
for ax, (s, sn) in zip(axs, SOL):
    x = [r[f"FC_{s}"] for r in R]; y = [r[f"SPECTRA_{s}"] for r in R]; c = [SZ[r["cluster"]] for r in R]
    sc = ax.scatter(x, y, c=c, cmap="viridis", norm=norm, s=6, lw=0, alpha=0.8, label="one task (colour = network size)")
    lo, hi = 5, 2e4; ax.plot([lo, hi], [lo, hi], color="#444", lw=0.8, label="same time"); ax.plot([lo, hi], [2 * lo, 2 * hi], ":", color="#888", lw=0.7, label="2x slower / 2x faster")
    ax.plot([lo, hi], [lo / 2, hi / 2], ":", color="#888", lw=0.7)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    nf = sum(b > a for a, b in zip(x, y)); ax.set_title(f"{sn}\nFasterCC faster in {nf:,} of {len(x):,} tasks")
    ax.set_xlabel("FasterCC time per task (s)"); ax.text(0.03, 0.95, "above the line: SPECTRA slower", transform=ax.transAxes, va="top", fontsize=7, color="#444")
    rows += [[r["task"], sn, r[f"FC_{s}"], r[f"SPECTRA_{s}"], SZ[r["cluster"]], r["level"]] for r in R]
axs[0].set_ylabel("SPECTRA time per task (s)")
for ax in axs: ax.legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(0, 0.9), markerscale=2)
cb = fig.colorbar(sc, ax=axs, pad=0.02, shrink=0.9); cb.set_label("network size (reactions)")
save(fig, "figSP1_task_time_fastercc_vs_spectra", rows, ["task", "solver", "FasterCC_s", "SPECTRA_s", "reactions", "blocked"], f"""
**Figure SP1 | Time per task, FasterCC against SPECTRA (1,140 cluster tasks per solver).** Each point is one task (one network
of 21,200 to 183,920 reactions with 0 to 500 reactions closed at random, one replicate); colour = network size. Solid line:
same time; dotted lines: 2x faster or slower. Points above the line are tasks where SPECTRA took longer than FasterCC. The 15
SPECTRA times from licence-failed runs are those of the SPEC15 reruns. Data: figSP1_task_time_fastercc_vs_spectra.csv; code:
make_spectra_figs.py.""")

# ---------------- SP2: heatmap of median SPECTRA/FasterCC ratio, cluster x blocked level
rows = []; fig, axs = plt.subplots(1, 2, figsize=(7.2, 4.2), sharey=True)
M = {}
for s, sn in SOL:
    M[s] = [[st.median(r[f"SPECTRA_{s}"] / r[f"FC_{s}"] for r in R if r["cluster"] == c and r["level"] == l) for l in LV] for c in CL]
vals = [v for s in M for row in M[s] for v in row]; nrm = TwoSlopeNorm(vcenter=0, vmin=min(math.log2(v) for v in vals), vmax=max(math.log2(v) for v in vals))
for ax, (s, sn) in zip(axs, SOL):
    im = ax.imshow([[math.log2(v) for v in row] for row in M[s]], cmap="PuOr_r", norm=nrm, aspect="auto")
    for i, c in enumerate(CL):
        for j, l in enumerate(LV):
            v = M[s][i][j]; ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.5, color="#111" if abs(math.log2(v)) < 1.2 else "#fff")
            rows.append([sn, c, SZ[c], l, round(v, 4)])
    ax.set_xticks(range(len(LV))); ax.set_xticklabels(LV); ax.set_xlabel("blocked reactions"); ax.set_title(sn); ax.grid(False)
axs[0].set_yticks(range(len(CL))); axs[0].set_yticklabels([f"c{c:02d} ({SZ[c]:,})" for c in CL]); axs[0].set_ylabel("cluster network (reactions)")
cb = fig.colorbar(im, ax=axs, pad=0.02, shrink=0.9); cb.set_label("colour: log2(SPECTRA time / FasterCC time)\norange (> 0): FasterCC faster; purple (< 0): SPECTRA faster")
fig.text(0.5, -0.02, "Number in each cell = median ratio SPECTRA time / FasterCC time over the 10 replicates (1.00 = same time).", ha="center", fontsize=7)
save(fig, "figSP2_spectra_ratio_heatmap", rows, ["solver", "cluster", "reactions", "blocked", "median_ratio_SPECTRA_over_FasterCC"], """
**Figure SP2 | Where FasterCC or SPECTRA is faster.** Median over 10 replicates of the per-task ratio SPECTRA time / FasterCC
time, for every network (rows) and number of blocked reactions (columns). Orange (ratio > 1): FasterCC faster; purple (ratio
< 1): SPECTRA faster; the number in each cell is the ratio. Data: figSP2_spectra_ratio_heatmap.csv; code: make_spectra_figs.py.""")

# ---------------- groups of real models
def group(m):
    if m.startswith("Fetus_"): return "Fetal tissues"
    if m.startswith(("Harvetta", "Harvey")): return "Whole-body"
    if m in ("Microbiota",): return "Microbiota"
    if m.startswith("ME_"): return "ME model"
    if m.startswith(("Recon", "Human")): return "Human generic"
    if m in ("DAS", "subDas", "cardiac_mit_glcuptake_atpmax", "ecoli_core_model"): return "Small test"
    if m.startswith(("Ec_", "STM", "Sc_", "iIT", "iJO", "iML", "iMM")): return "Curated microbial"
    return "AGORA microbes"
GORD = ["Small test", "AGORA microbes", "Curated microbial", "Human generic", "Fetal tissues", "Whole-body", "Microbiota", "ME model"]
f = lambda x: float(x) if x not in ("", "NaN", "nan", "NA") else float("nan")
for m in RM: m["group"] = group(m["model"]); m["in_total"] = int(m["in_total"])
RMs = sorted(RM, key=lambda m: (GORD.index(m["group"]), m["in_total"]))
TN = {"FasterCC": "NEW", "FASTCC": "FASTCC", "SPECTRA": "SPECTRA"}

# ---------------- SP3: real-model run times
rows = []; fig, ax = plt.subplots(figsize=(6.6, 9.0)); y = list(range(len(RMs)))[::-1]
for meth in ("FasterCC", "FASTCC", "SPECTRA"):
    t = [f(m[f"{TN[meth]}_time_s"]) for m in RMs]
    ax.scatter([v for v in t if v == v], [yy for yy, v in zip(y, t) if v == v], marker=MK[meth], s=16, color=COL[meth], label=meth, zorder=3, lw=0)
    rows += [[m["model"], m["group"], m["in_total"], meth, v] for m, v in zip(RMs, t)]
ax.set_yticks(y); ax.set_yticklabels([f"{m['model'][:30]} ({m['in_total']:,})" for m in RMs], fontsize=5.5); ax.set_xscale("log")
ax.set_xlabel("run time (s, log scale)"); ax.set_ylim(-1, len(RMs)); ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.03), ncol=3, frameon=False)
prev = None
for yy, m in zip(y, RMs):
    if m["group"] != prev:
        ax.axhline(yy + 0.5, color="#bbb", lw=0.5); ax.text(1.0, yy + 0.25, m["group"], transform=ax.get_yaxis_transform(), ha="right", va="top", fontsize=6, style="italic", color="#555")
    prev = m["group"]
gtab = []
for g in GORD:
    ms = [m for m in RM if m["group"] == g]; r3 = [m for m in ms if all(f(m[f"{TN[k]}_time_s"]) > 0 for k in TN)]
    if g == "Whole-body":
        r2 = [m for m in r3 if "aafix" in m["model"]]
        gtab.append(["Whole-body, aafix only (feasible)", len(r2)] + [round(gm([f(m[f"{TN[k]}_time_s"]) for m in r2]), 2) for k in TN] +
                    [sum(f(m["SPECTRA_time_s"]) > f(m["NEW_time_s"]) for m in r2), sum(f(m["FASTCC_time_s"]) > f(m["NEW_time_s"]) for m in r2)])
    if r3: gtab.append([g, len(r3)] + [round(gm([f(m[f"{TN[k]}_time_s"]) for m in r3]), 2) for k in TN] +
                       [sum(f(m["SPECTRA_time_s"]) > f(m["NEW_time_s"]) for m in r3), sum(f(m["FASTCC_time_s"]) > f(m["NEW_time_s"]) for m in r3)])
save(fig, "figSP3_real_model_times", rows, ["model", "group", "reactions", "method", "time_s"], """
**Figure SP3 | Run time of FasterCC, FASTCC and SPECTRA on the 58 real models.** Models grouped (italic labels) and ordered by
size within each group; reactions in brackets. FasterCC has no time for ME_dbl_E0_L1_00e_2 (its run is the ME restart still
running). Times as stored in results_real (one run per model and method). On small and mid-size real models SPECTRA is
usually the fastest; FasterCC is faster than FASTCC on the larger models. The FasterCC times of Harvey_1_03c and Harvetta_1_03c (about 23 h) are dominated by careful-mode re-solves: these two models are infeasible with their own bounds at 1e-9, so about 850 LPs fall back to single-thread primal simplex (tableR6b_careful_mode_1_03c); their feasible aafix versions take 47-54 min. Per-group geometric means: tableSP2_real_times_by_group.
Data: figSP3_real_model_times.csv; code: make_spectra_figs.py.""")
table("tableSP2_real_times_by_group", ["group", "models (all 3 finished)", "FasterCC geo-mean s", "FASTCC geo-mean s", "SPECTRA geo-mean s",
                                       "models SPECTRA slower than FasterCC", "models FASTCC slower than FasterCC"], gtab, """
**Table SP2 | Real models: run time per group (geometric mean over the models where all three methods finished).** Source:
methods_real_16Nov_models.csv (results_real).""")

# ---------------- SP4: agreement with FVA on real models
rows = []; err = [m for m in RMs if any(f(m.get(f"{TN[k]}_{e}", "")) > 0 for k in TN for e in ("falseBlocked", "falseConsistent"))]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 0.17 * len(err) + 1.2), sharey=True); yy = list(range(len(err)))[::-1]
for ax, (e, title) in zip(axs, (("falseBlocked", "FVA-consistent reactions removed"), ("falseConsistent", "kept although FVA calls them blocked"))):
    for j, meth in enumerate(("FasterCC", "FASTCC", "SPECTRA")):
        v = [f(m[f"{TN[meth]}_{e}"]) for m in err]
        ax.barh([a + (1 - j) * 0.27 for a in yy], [x if x == x else 0 for x in v], height=0.26, color=COL[meth], label=meth)
        rows += [[m["model"], meth, e, x] for m, x in zip(err, v)]
    ax.set_xscale("symlog", linthresh=1); ax.set_title(title); ax.set_xlabel("reactions (symmetric log)")
axs[0].set_yticks(yy); axs[0].set_yticklabels([m["model"][:30] for m in err], fontsize=6)
for ax in axs: ax.legend(frameon=False, loc="upper right", fontsize=7)
save(fig, "figSP4_real_models_vs_fva", rows, ["model", "method", "error_type", "reactions"], f"""
**Figure SP4 | Real models: where each method's kept set differs from FVA.** Only the {len(err)} models with any difference are
shown (the others agree with FVA for all three methods). Left: reactions FVA calls consistent that the method removed (an error
that loses network function). Right: reactions FVA calls blocked that the method kept. FasterCC's kept-but-blocked reactions in
the fetal, curated microbial and AGORA models are the guard-kept reactions; ME counts are against a floating-point FVA (the
exact checks of 14-15 Nov are the reference for ME). FVA as stored by the benchmark in each result file. Data:
figSP4_real_models_vs_fva.csv; code: make_spectra_figs.py.""")

# ---------------- SP5: reversible reactions left in the output (downstream difficulty)
rows = []; sel = [m for m in RMs if f(m["NEW_n"]) > 0 and int(m["in_rev"]) > 0]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 7.6), gridspec_kw={"width_ratios": [1.25, 1]}); ax = axs[0]; yy = list(range(len(sel)))[::-1]
for meth, key, mk in (("SPECTRA", "SPECTRA_rev", "^"), ("FASTCC", "FASTCC_rev", "s"), ("FasterCC", "NEWout_rev", "o")):
    p = [100 * f(m[key]) / f(m[f"{TN[meth]}_n"]) if f(m[f"{TN[meth]}_n"]) > 0 else float("nan") for m in sel]
    ax.scatter(p, yy, marker=mk, s=34 if meth == "SPECTRA" else 12, facecolors="none" if meth == "SPECTRA" else COL[meth], edgecolors=COL[meth], lw=0.8 if meth == "SPECTRA" else 0, zorder=3 if meth == "SPECTRA" else 4, label=meth + (" (after its reorientation)" if meth == "FasterCC" else " (input directions)"))
    rows += [[m["model"], meth, int(f(m[key])), int(f(m[f"{TN[meth]}_n"])), round(x, 2)] for m, x in zip(sel, p)]
ax.set_yticks(yy); ax.set_yticklabels([m["model"][:28] for m in sel], fontsize=5.5); ax.set_xlabel("reversible reactions, % of reactions kept")
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.07), ncol=1, fontsize=6.5); ax.set_title("A  real models", loc="left"); ax.set_xlim(0, 65)
ax = axs[1]; bt = collections.defaultdict(list)
for t in BT:
    if t["folder"] == "results_full_gurobi_spectra" and t["n"] not in ("", "0"):
        bt[(int(t["cluster"]), t["method"])].append(t)
for meth, key, mk in (("SPECTRA", "kept_rev", "^"), ("FASTCC", "kept_rev", "s"), ("FasterCC", "nRev_reported", "o")):
    tm = {"FasterCC": "NEW"}.get(meth, meth)
    v = [st.median(100 * float(t[key]) / float(t["n"]) for t in bt[(c, tm)]) for c in CL]
    ax.plot([SZ[c] for c in CL], v, marker=mk, ms=3.5, lw=1, color=COL[meth], label=meth)
    rows += [[f"cluster c{c:02d}", meth + " (median, Gurobi tasks)", "", "", round(x, 2)] for c, x in zip(CL, v)]
ax.set_xlabel("cluster network size (thousand reactions)"); ax.set_ylim(0, 55); ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v/1000:.0f}")); ax.set_ylabel("reversible reactions, % of kept (median of 60 tasks)")
ax.set_title("B  cluster benchmark (Gurobi)", loc="left"); ax.legend(frameon=False, fontsize=6.5, loc="lower left")
save(fig, "figSP5_reversible_in_output", rows, ["model_or_cluster", "method", "reversible_kept", "kept", "percent_reversible"], """
**Figure SP5 | Reversible reactions left in the consistent network.** FASTCC and SPECTRA return the kept reactions with the
input bounds, so every reaction reversible in the input stays reversible. FasterCC also returns the bounds after its
reorientation (fixIrr turns backward-only reactions forward; the structure analysis makes reversible reactions irreversible
where the network and constraints allow only one direction), so fewer reversible reactions remain for downstream analysis
(sampling, FVA, loop removal, MILP methods, where each reversible reaction adds work). A: % of kept reactions that are
reversible, per real model (reversible = lb < 0 < ub; SPECTRA, open triangles, keeps the same reversible share as FASTCC and is mostly hidden under it; the same in B). B: same, median over the 60 Gurobi tasks of each cluster network
(FasterCC: res.nRev after reorientation). FasterCC's saved output model keeps the input orientation and bounds; the
reoriented counts come from its returned bounds. Data: figSP5_reversible_in_output.csv; code: make_spectra_figs.py.""")

# ---------------- Table SP1: benchmark agreement incl. SPEC15
n = collections.defaultdict(dict)
for t in BT:
    key = (int(t["cluster"]), int(t["level"]), int(t["rep"]))
    if t["folder"] == "results_full_cplex" and t["method"] == "NEW": n[key][("CPLEX", "FasterCC")] = int(t["n"] or 0)
    if t["folder"] == "results_full_cplex_spectra" and t["method"] in ("FASTCC", "SPECTRA"): n[key][("CPLEX", t["method"])] = int(t["n"] or 0)
    if t["folder"] == "results_full_gurobi_spectra": n[key][("Gurobi", {"NEW": "FasterCC"}.get(t["method"], t["method"]))] = int(t["n"] or 0)
rerun = 0
for t in S15:
    key = (int(t["cluster"]), int(t["level"]), int(t["rep"])); sol = "CPLEX" if t["solver"] == "cplex" else "Gurobi"
    if t["method"] == "SPECTRA" and n[key].get((sol, "SPECTRA"), 0) == 0: n[key][(sol, "SPECTRA")] = int(t["n"]); rerun += 1
    if t["method"] == "FASTCC" and n[key].get((sol, "FASTCC"), 0) == 0: n[key][(sol, "FASTCC")] = int(t["n"])
trows = []
for sol in ("CPLEX", "Gurobi"):
    for a, b in (("FasterCC", "FASTCC"), ("FasterCC", "SPECTRA"), ("FASTCC", "SPECTRA")):
        pairs = [(d.get((sol, a), 0), d.get((sol, b), 0)) for d in n.values()]; both = [(x, y) for x, y in pairs if x > 0 and y > 0]
        trows.append([sol, f"{a} vs {b}", len(both), sum(x == y for x, y in both), sum(x > y for x, y in both), sum(x < y for x, y in both), len(pairs) - len(both)])
table("tableSP1_kept_count_agreement", ["solver", "comparison", "tasks with both results", "same number kept", "first keeps more", "second keeps more", "missing"], trows, f"""
**Table SP1 | Cluster benchmark: number of reactions kept, method against method, same task.** After the {rerun} SPEC15 SPECTRA
reruns replace the licence-failed results (n = 0). Differences of one reaction are biomass reactions (see
methods_bench_16Nov_diff_ids_finished.csv); FasterCC CPLEX +1 in task 1119 is the KEPTCHK case of 13 Nov. Source:
methods_bench_tasks.csv + spec15_tasks.tsv.""")
