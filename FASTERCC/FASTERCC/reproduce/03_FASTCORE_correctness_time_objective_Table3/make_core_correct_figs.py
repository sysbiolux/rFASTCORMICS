#!/usr/bin/env python3
"""make_core_correct_figs.py -- preprint Figure 4, Table 2, Supplementary Table 2 and Table 3; ledger Figure R3, Table R7
(FASTCORE with and without FasterCC preprocessing), Figure R5 (correctness), Table R6 (elapsed and CPU time on the largest
models) and Table R22 (objective in each method's result), rebuilt from the raw benchmark outputs.
Inputs (data/):
  CORE_by_task.csv            250 FASTCORE tasks (compareA/compare_A.m on results_full_cplex/CORE; clusters 1-10, CPLEX)
  methods_real_16Nov_models.csv, methods_bench_tasks.csv   kept sets / FVA errors (export_compare.m + compare_methods.py)
  fva_check/*.json            independent FVA consistency check of the cluster networks (fcc_consistency_check_v2.py,
                              Gurobi 12; one file per network; counts only)
  sacct_largest.txt, tasks_real.tsv   Slurm accounting of the largest real-model runs; task index -> model, method
  BIGOBJ_16226422.out, BIGOBJ2_16210166.out   objective tests on the largest networks
  table3/                     FVA of iMM1415, iMM1865 (xlsx), Recon3D_mat, HumanGEM12971 (csv); FasterCC 6 Nov pipeline
                              runs of the two human models: per-step counts and every reoriented reaction
Writes into OUT (refuses to overwrite). Usage: python3 make_core_correct_figs.py data OUT"""
import csv, glob, json, math, os, re, statistics as st, sys, collections
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
P = lambda *a: os.path.join(src, *a)
COL = {"FasterCC": "#0072B2", "FASTCC": "#E69F00", "SPECTRA": "#009E73", "raw": "#E69F00", "clean": "#0072B2", "grey": "#8c8c8c"}
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
q = lambda v, k: st.quantiles(v, n=4)[k]
gm = lambda v: math.exp(sum(math.log(x) for x in v) / len(v))

# ======================= FASTCORE
C = list(csv.DictReader(open(P("CORE_by_task.csv"))))
for r in C:
    for k in ("cluster", "core", "rep", "nRxns"): r[k] = int(r[k])
    for k in ("time_preclean", "time_noPre", "time_pre", "F1exact_noPre", "F1exact_pre", "jaccard"): r[k] = float(r[k])
    r["speedup"] = r["time_noPre"] / r["time_pre"]
CORES = sorted({r["core"] for r in C}); CLU = sorted({r["cluster"] for r in C}); SZ = {r["cluster"]: r["nRxns"] for r in C}
assert len(C) == 250 and all(r["coreKept_noPre"] == "1" and r["coreKept_pre"] == "1" for r in C)
rows = []; fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
ax = axs[0]
for key, lab, colr, ls in (("time_noPre", "without preprocessing", COL["raw"], "--"), ("time_pre", "with FasterCC preprocessing", COL["clean"], "-")):
    med = [st.median(r[key] for r in C if r["core"] == c) for c in CORES]; lo = [q([r[key] for r in C if r["core"] == c], 0) for c in CORES]
    hi = [q([r[key] for r in C if r["core"] == c], 2) for c in CORES]
    ax.fill_between(range(len(CORES)), lo, hi, color=colr, alpha=0.18, lw=0); ax.plot(range(len(CORES)), med, ls, color=colr, marker="o", ms=3, lw=1.4, label=lab)
    rows += [["time", lab, c, round(m, 2), round(a, 2), round(b, 2)] for c, m, a, b in zip(CORES, med, lo, hi)]
med = [st.median(r["time_pre"] + r["time_preclean"] for r in C if r["core"] == c) for c in CORES]
ax.plot(range(len(CORES)), med, ":", color=COL["grey"], marker="s", ms=2.5, lw=1, label="cleaning + FASTCORE (one core set)")
rows += [["time", "cleaning + FASTCORE", c, round(m, 2), "", ""] for c, m in zip(CORES, med)]
ax.set_xticks(range(len(CORES))); ax.set_xticklabels(CORES); ax.set_xlabel("core reactions"); ax.set_ylabel("FASTCORE time (s), median and IQR")
ax.legend(frameon=False, fontsize=6.5); ax.set_title("a", loc="left", fontweight="bold")
ax = axs[1]
med = [st.median(r["speedup"] for r in C if r["core"] == c) for c in CORES]; lo = [q([r["speedup"] for r in C if r["core"] == c], 0) for c in CORES]
hi = [q([r["speedup"] for r in C if r["core"] == c], 2) for c in CORES]
ax.fill_between(range(len(CORES)), lo, hi, color=COL["clean"], alpha=0.18, lw=0); ax.plot(range(len(CORES)), med, "-", color=COL["clean"], marker="o", ms=3, lw=1.4, label="median per-task speed-up (band: IQR)")
ax.axhline(1, ls=":", color="#444", lw=0.8, label="no gain (same time)"); ax.legend(frameon=False, fontsize=6.5, loc="lower right"); ax.set_xticks(range(len(CORES))); ax.set_xticklabels(CORES); ax.set_xlabel("core reactions")
ax.set_ylabel("speed-up (time without / with preprocessing)"); ax.set_title("b", loc="left", fontweight="bold"); ax.set_ylim(0, None)
rows += [["speedup", "per task", c, round(m, 3), round(a, 3), round(b, 3)] for c, m, a, b in zip(CORES, med, lo, hi)]
n10 = sum(r["speedup"] >= 1 / 0.9 for r in C); nfast = sum(r["speedup"] > 1 for r in C)
save(fig, "figP4_R3_fastcore_preprocessing", rows, ["panel", "series", "core_reactions", "median", "q1", "q3"], f"""
**Figure 4 (preprint) / R3 (ledger) | FASTCORE with and without FasterCC preprocessing.** a, Median (line) and interquartile
range (band) of FASTCORE runtime on the raw expanded cluster network and on the same network after FasterCC cleaning, by
core size (50 tasks per core size: clusters 1-10 of 12,160 to 93,520 reactions x 5 replicates; CPLEX); dotted grey, cleaning
plus one FASTCORE run. b, Per-task speed-up (time without / with preprocessing), median and IQR; dotted line, no gain.
FASTCORE was faster on the cleaned network in {nfast} of 250 tasks and at least 10% faster (time ratio >= 1/0.9) in {n10}.
The core was kept in all 250 tasks with both inputs; the two outputs were never identical (Jaccard index in
tableR7_fastcore_by_cluster). Note: the preprint's networks started at 10,600 reactions; these expanded inputs start at
12,160 (cluster 1). Data: figP4_R3_fastcore_preprocessing.csv; source CORE_by_task.csv; code: make_core_correct_figs.py.""")
# Supplementary Table 2: mean of 5 replicates per cluster x core size
srows = []
for c in CLU:
    for k in CORES:
        t = [r for r in C if r["cluster"] == c and r["core"] == k]
        a, b = st.mean(r["time_noPre"] for r in t), st.mean(r["time_pre"] for r in t)
        srows.append([SZ[c], k, round(a, 1), round(st.stdev(r["time_noPre"] for r in t), 1), round(b, 1), round(st.stdev(r["time_pre"] for r in t), 1),
                      round(a / b, 2), sum(r["speedup"] > 1 for r in t), round(st.mean(r["time_preclean"] for r in t), 1),
                      round(st.mean(r["F1exact_noPre"] for r in t), 3), round(st.mean(r["F1exact_pre"] for r in t), 3)])
table("tableS2_fastcore_by_size_and_core", ["input reactions", "core reactions", "FASTCORE mean s (no preprocessing)", "SD", "FASTCORE mean s (FasterCC preprocessing)", "SD",
                                            "speed-up (ratio of means)", "replicates faster (of 5)", "FasterCC cleaning mean s", "F1 exact, raw", "F1 exact, cleaned"], srows, """
**Supplementary Table 2 | FASTCORE running time with and without FasterCC preprocessing, every network size and core size.**
Mean (SD) over 5 replicates; speed-up = mean time without / mean time with preprocessing; the cleaning time is paid once per
network. Source: CORE_by_task.csv (CPLEX).""")
top = sorted(srows, key=lambda r: -r[6])[:5]
table("table2_fastcore_top5", ["input reactions", "core reactions", "time, no preprocessing (s)", "time with FasterCC preprocessing (s)", "speed-up"],
      [[r[0], r[1], r[2], r[4], r[6]] for r in top], """
**Table 2 | FASTCORE running time improvement due to preprocessing with FasterCC: the five largest speed-ups** (ratio of mean
times over 5 replicates; full table: tableS2_fastcore_by_size_and_core). Source: CORE_by_task.csv.""")
rrows = []
for c in CLU:
    t = [r for r in C if r["cluster"] == c]
    rrows.append([f"c{c:02d}", SZ[c], round(st.median(r["time_preclean"] for r in t), 1), round(st.median(r["time_noPre"] for r in t), 1),
                  round(st.median(r["time_pre"] for r in t), 1), round(st.median(r["speedup"] for r in t), 2),
                  f"{st.median(r['F1exact_noPre'] for r in t):.3f} -> {st.median(r['F1exact_pre'] for r in t):.3f}", round(st.median(r["jaccard"] for r in t), 3)])
table("tableR7_fastcore_by_cluster", ["cluster", "reactions", "FasterCC cleaning s", "FASTCORE raw s", "FASTCORE cleaned s", "speed-up (median per task)", "F1 exact raw -> cleaned", "Jaccard raw vs cleaned"], rrows, """
**Table R7 | FASTCORE with and without FasterCC preprocessing, by cluster.** Medians of 25 tasks per cluster (5 core sizes x 5
replicates). Source: CORE_by_task.csv.""")

# ======================= R5: correctness
RM = list(csv.DictReader(open(P("methods_real_16Nov_models.csv"))))
f = lambda x: float(x) if x not in ("", "NaN", "nan", "NA") else float("nan")
rows = []; fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2))
ax = axs[0]; pts = [m for m in RM if f(m["NEW_time_s"]) > 0 and f(m["SPECTRA_time_s"]) > 0]
for m in pts:
    k = f(m["SPECTRA_falseConsistent"]); ax.scatter(f(m["NEW_time_s"]), f(m["SPECTRA_time_s"]), s=8 + 3 * math.sqrt(k), color=COL["SPECTRA"], alpha=0.7, lw=0)
    rows.append(["a", m["model"], f(m["NEW_time_s"]), f(m["SPECTRA_time_s"]), int(k)])
    off = {"Microbiota": (-10, 8), "Harvetta_aafix_E6": (-60, 10), "Harvey_aafix_E6": (-58, -12), "Harvey_1_03c": (-40, -14), "Harvetta_1_03c": (-46, 8)}
    if k >= 100: ax.annotate(f"{m['model']} ({int(k)})", (f(m["NEW_time_s"]), f(m["SPECTRA_time_s"])), fontsize=5.5, xytext=off.get(m["model"], (3, -8)), textcoords="offset points",
                             arrowprops=dict(arrowstyle="-", lw=0.4, color="#666"))
ax.set_xscale("log"); ax.set_yscale("log"); lim = (0.5, 2e5); ax.plot(lim, lim, "--", color="#444", lw=0.8, label="same time (below: SPECTRA faster)")
for kk in (0, 100, 1000): ax.scatter([], [], s=8 + 3 * math.sqrt(kk), color=COL["SPECTRA"], alpha=0.7, lw=0, label=f"one model; SPECTRA keeps {kk} FVA-blocked")
ax.legend(frameon=False, fontsize=5.5, loc="upper left"); ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("FasterCC time (s)"); ax.set_ylabel("SPECTRA time (s)"); ax.set_title(f"a  real models ({len(pts)} with both results)", loc="left")
ax = axs[1]
J = {}
for fn in glob.glob(P("fva_check", "*.json")):
    d = json.load(open(fn)); m = re.search(r"c(\d+)_l(\d+)_r(\d+)", d["tag"]); J[tuple(int(x) for x in m.groups())] = d
N = collections.defaultdict(dict)
for t in csv.DictReader(open(P("methods_bench_tasks.csv"))):
    key = (int(t["cluster"]), int(t["level"]), int(t["rep"]))
    if t["n"] in ("", "0"): continue
    if t["folder"] == "results_full_cplex" and t["method"] == "NEW": N[key][("CPLEX", "FasterCC")] = int(t["n"])
    if t["folder"] == "results_full_cplex_spectra" and t["method"] in ("FASTCC", "SPECTRA"): N[key][("CPLEX", t["method"])] = int(t["n"])
    if t["folder"] == "results_full_gurobi_spectra": N[key][("Gurobi", {"NEW": "FasterCC"}.get(t["method"], t["method"]))] = int(t["n"])
both = sorted(k for k in J if k in N); res = collections.defaultdict(lambda: [0, 0, 0, 0])   # tasks, kept more, kept fewer, extra total
for k in both:
    ref = J[k]["n_consistent"]
    for (s, meth), n in N[k].items():
        x = res[(s, meth)]; x[0] += 1; x[1] += n > ref; x[2] += n < ref; x[3] += max(n - ref, 0)
meths = ["FasterCC", "FASTCC", "SPECTRA"]; w = 0.38
for j, s in enumerate(("CPLEX", "Gurobi")):
    v = [res[(s, m)][1] for m in meths]
    ax.bar([i + (j - 0.5) * w for i in range(3)], v, width=w, color=[COL[m] for m in meths], hatch=None if s == "CPLEX" else "///", edgecolor="white", lw=0.5,
           label=s)
    for i, m in enumerate(meths): rows.append(["b", f"{s} {m}", res[(s, m)][0], res[(s, m)][1], res[(s, m)][2], res[(s, m)][3]])
ax.set_xticks(range(3)); ax.set_xticklabels(meths); ax.set_ylabel(f"tasks keeping more reactions than FVA\n(of {len(both)} with an FVA reference)")
ax.set_title("b  cluster tasks vs independent FVA", loc="left")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor="#999", label="CPLEX"), Patch(facecolor="#999", hatch="///", edgecolor="white", label="Gurobi")], frameon=False, fontsize=6.5)
fewer = sum(res[k][2] for k in res)
save(fig, "figR5_correctness", rows, ["panel", "model_or_series", "x_or_tasks", "y_or_kept_more", "kept_blocked_or_kept_fewer", "extra_reactions_total"], f"""
**Figure R5 | Correctness of FasterCC, SPECTRA and FASTCC.** a, Real models: SPECTRA time against FasterCC time (log-log;
dashed, equal time; the two points near 8e4 s on the x axis are Harvey_1_03c and Harvetta_1_03c, slowed by careful-mode re-solves, see tableR6b_careful_mode_1_03c); marker area grows with the number of FVA-blocked reactions SPECTRA keeps (labelled when >= 100). b,
Cluster tasks with an independent FVA reference ({len(both)} tasks: the networks of clusters 2-6 that are in both the speed
benchmark and the FVA consistency check, which so far covers clusters 1-6 and two networks of cluster 7; Gurobi 12, threshold
0.99 eps, every LP optimal): number of tasks where a method keeps more reactions than FVA finds consistent; solid CPLEX,
hatched Gurobi. Tasks where a method keeps fewer than FVA (removing a reaction that can carry flux): {fewer} in total. This
comparison uses counts (the FVA check stores counts, not reaction IDs); the reaction IDs of the differences between methods are
in methods_bench_16Nov_diff_ids_finished.csv (all biomass reactions). Data: figR5_correctness.csv; code: make_core_correct_figs.py.""")

# ======================= Table R6: elapsed and CPU time (Slurm accounting)
T = {r["task"]: r for r in csv.DictReader(open(P("tasks_real.tsv")), delimiter="\t")}
def sec(s):
    d = 0
    if "-" in s: d, s = s.split("-"); d = int(d)
    p = [float(x) for x in s.split(":")]
    while len(p) < 3: p = [0] + p
    return d * 86400 + p[0] * 3600 + p[1] * 60 + p[2]
A = collections.defaultdict(dict)
for r in csv.DictReader(open(P("sacct_largest.txt")), delimiter="|"):
    if "." in r["JobID"]: continue
    t = T[r["JobID"].split("_")[1]]; A[t["model"]][{"NEW": "FasterCC"}.get(t["method"], t["method"])] = (sec(r["Elapsed"]), sec(r["TotalCPU"]), int(r["AllocCPUS"]), r["JobID"], r["State"])
trows = []
for m in sorted(A, key=lambda m: -A[m]["FASTCC"][0]):
    a = A[m]; fc = a["FasterCC"]
    trows.append([m, fc[2], *[round(a[k][0]) for k in ("FasterCC", "FASTCC", "SPECTRA")], *[round(a[k][1]) for k in ("FasterCC", "FASTCC", "SPECTRA")],
                  " / ".join(f"{a[k][1] / a[k][0]:.1f}" for k in ("FasterCC", "FASTCC", "SPECTRA")),
                  f"{a['FASTCC'][0] / fc[0]:.1f} -> {a['FASTCC'][1] / fc[1]:.1f}", f"{a['SPECTRA'][0] / fc[0]:.2f} -> {a['SPECTRA'][1] / fc[1]:.2f}",
                  " ".join(a[k][3] for k in ("FasterCC", "FASTCC", "SPECTRA"))])
table("tableR6_wall_and_cpu_time", ["model", "CPUs allocated", "elapsed s FasterCC", "elapsed s FASTCC", "elapsed s SPECTRA", "CPU s FasterCC", "CPU s FASTCC", "CPU s SPECTRA",
                                    "cores used FasterCC / FASTCC / SPECTRA", "FASTCC / FasterCC: elapsed -> CPU", "SPECTRA / FasterCC: elapsed -> CPU", "Slurm jobs"], trows, """
**Table R6 | Elapsed and CPU time on the largest real models (Slurm accounting).** Elapsed and TotalCPU of each run's Slurm
job (sacct), cores used = CPU / elapsed. Each method ran in its own job with the same allocation per model class (38 CPUs
for the whole-body, Microbiota and ME models, 10 for the human models). The ME FasterCC run completed in Slurm but has no
exported result (ME restart), so its time is listed but not compared elsewhere. Sources: sacct_largest.txt, tasks_real.tsv.""")

# ======================= Table R22: objective
def parse(fn):
    D = collections.OrderedDict(); cur = None
    for line in open(fn):
        m = re.match(r"#+ (\S+): (\d+) reactions.*objective (\S+) \((?:model )?bounds ([^)]*)\)", line)
        if m: cur = m.group(1); D[cur] = {"obj": m.group(3), "bounds": m.group(4), "rx": int(m.group(2))}; continue
        m = re.match(r"== (\w+)\s+kept\s+(\d+) \| objective kept (\d) \| FEASIBLE with the model bounds: (\w+) \| largest objective flux \(bounds opened\): (\S+)", line)
        if m and cur: D[cur][m.group(1)] = (int(m.group(2)), m.group(3), m.group(4), m.group(5))
    return D
D = parse(P("BIGOBJ_16226422.out")); orows = []
for net, d in D.items():
    orows.append([net, f"{d['obj']} ({d['bounds']})"] + [f"{d[k][0]:,}; feasible {d[k][2]}; max {d[k][3]}" if k in d else "no result" for k in ("FULL", "NEW", "FASTCC", "SPECTRA")])
table("tableR22_objective", ["network", "objective (bounds)", "full input network", "FasterCC", "FASTCC", "SPECTRA"], orows, """
**Table R22 | Objective in each method's result on the largest networks.** For each network the objective reaction is tested
in the full input and in each method's result: reactions kept; feasible with the model's bounds; largest objective flux with
the bounds opened. Job BIGOBJ 16226422 (11 Nov, all four networks with results for all methods). The whole-body _1_03c
models are infeasible with their own bounds already in the full network (BIGOBJ2 16210166: feasible only at feasibility
tolerance 1e-6), so 'max 0' there reflects the model, not the method; on Microbiota only FasterCC's network can carry
community biomass (459.2, as the full network). Source: BIGOBJ_16226422.out, BIGOBJ2_16210166.out.""")

# ======================= Table 3: reversible reactions by bounds, FasterCC preprocessing, FVA
import openpyxl
t3 = []; eps = 1e-6
for name, fn in (("iMM1415", "FVA_iMM1415_fullmodel.xlsx"), ("iMM1865", "FVA_iMM1865_fullmodel.xlsx")):
    ws = openpyxl.load_workbook(P("table3", fn), read_only=True)["full_model"]; R = list(ws.iter_rows(min_row=2, values_only=True))
    t3.append([name, len(R), sum(1 for r in R if r[4] < 0 < r[5]), "not run with the current version", "", sum(1 for r in R if r[7] < -eps and r[8] > eps)])
ST = {r["model"]: r for r in csv.DictReader(open(P("table3", "steps_Recon3D_mat_HumanGEM12971.csv")))}
RO = list(csv.DictReader(open(P("table3", "reoriented_Recon3D_mat_HumanGEM12971.csv"))))
chk = []
for name, label in (("Recon3D_mat", "Recon3D (10,600)"), ("HumanGEM12971", "Human-GEM (12,971)")):
    F = {r["rxn_id"]: (float(r["lb"]), float(r["ub"])) for r in csv.DictReader(open(P("table3", f"{name}_FVA_lbub.csv")))}
    ro = [r for r in RO if r["model"] == name and r["mode"] == "def"]; miss = sum(r["rxn"] not in F for r in ro)
    toirr = [r for r in ro if r["input_class"] == "rev" and r["lp_class"] != "rev"]
    agree = sum((F[r["rxn"]][1] <= eps) if r["flipped"] == "1" else (F[r["rxn"]][0] >= -eps) for r in toirr)
    fva_bi = sum(1 for v in F.values() if v[0] < -eps and v[1] > eps)
    t3.append([label, len(F), int(ST[name]["in_rev"]), int(ST[name]["lp_rev"]), f"{agree} of {len(toirr)}", fva_bi])
    chk.append([label, len(ro), miss, len(toirr), agree, len(toirr) - agree, int(ST[name]["lp_rev"]) - fva_bi])
table("table3_reversible_reactions", ["model", "reactions", "reversible (bounds)", "reversible after FasterCC preprocessing (no optional check)",
                                      "preprocessing reorientations confirmed by FVA", "reversible (FVA)"], t3, """
**Table 3 | Number of reversible reactions by bounds, after FasterCC preprocessing, and by FVA.** Reversible by bounds: lb < 0 <
ub. Reversible by FVA: minimum flux < -1e-6 and maximum flux > 1e-6 (the same counts at 1e-9 and 1e-4). FasterCC column: the
6 Nov pipeline runs (default mode) WITHOUT FasterCC's optional reversibility check, which was not run with the current version;
every reaction its preprocessing turned irreversible is irreversible in the same direction by FVA (column 5), and the
remaining difference to FVA is what the optional check is meant to close. The preprint's FasterCC column (equal to FVA) must
be reproduced with the optional check (planned after the 16 Nov maintenance, also on iMM1415 and iMM1865). The bounds and FVA
columns reproduce the preprint exactly. Sources: data/table3/.""")
table("table3_check_by_reaction_id", ["model", "reactions reoriented by FasterCC (any kind)", "IDs not found in the FVA file", "reversible -> irreversible",
                                      "same direction in FVA", "other direction in FVA", "left reversible by preprocessing minus reversible in FVA"], chk, """
**Table 3 check | FasterCC preprocessing against FVA, reaction by reaction.** Every reaction FasterCC reoriented is found by ID
in the FVA file of the same model (same reaction count), and every reversible reaction it made irreversible can carry flux
only in that direction by FVA. Source: data/table3/.""")
