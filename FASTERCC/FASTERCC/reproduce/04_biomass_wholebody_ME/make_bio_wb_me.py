#!/usr/bin/env python3
"""make_bio_wb_me.py -- ledger Table R23 (biomass reactions kept or removed wrongly), the SPECTRA whole-body table and the
ME snapshot table, with figures. Inputs (data/): biocheck/bio_cNN_fill2.csv (biomass check of 10 Nov: every biomass
reaction of every task decided by exact LPs, BIOCHK + fill-in, md5 of bio_fill.m 260ec994 in the ledger),
spectrawb_v1_table.txt (SPECTRA on the aafix whole-body models against FVA), me_fva_vs_quad.tsv (ME model: FasterCC's
state after the singleton stage (snapshot 4) per reaction, the FVA call and SPECTRA's result). Writes into OUT (no
overwrite). Usage: python3 make_bio_wb_me.py data OUT"""
import csv, glob, os, sys, collections
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
COL = {"FasterCC": "#0072B2", "FASTCC": "#E69F00", "SPECTRA": "#009E73", "FasterCC (preprint)": "#56B4E9"}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6",
                     "axes.axisbelow": True, "svg.fonttype": "none", "pdf.fonttype": 42})
def guard(*n):
    for x in n:
        if os.path.exists(os.path.join(out, x)): raise SystemExit(x + " exists")
def table(name, header, rows, cap):
    guard(name + ".csv", name + ".md")
    with open(os.path.join(out, name + ".csv"), "w", newline="") as h: w = csv.writer(h); w.writerow(header); w.writerows(rows)
    open(os.path.join(out, name + ".md"), "w").write(cap.strip() + "\n\n| " + " | ".join(header) + " |\n|" + "---|" * len(header) + "\n" +
        "\n".join("| " + " | ".join(f"{x:,}" if isinstance(x, int) and abs(x) > 999 else str(x) for x in r) + " |" for r in rows) + "\n"); print("written", name)
def save(fig, name, cap):
    guard(*[f"{name}.{e}" for e in ("svg", "pdf", "png", "md")])
    for e in ("svg", "pdf", "png"): fig.savefig(os.path.join(out, f"{name}.{e}"), dpi=300, bbox_inches="tight")
    plt.close(fig); open(os.path.join(out, name + ".md"), "w").write(cap.strip() + "\n"); print("written", name)

# ---- Biomass (Table R23 + figure)
R = []
for f in sorted(glob.glob(os.path.join(src, "biocheck", "bio_c*_fill2.csv"))): R += list(csv.DictReader(open(f)))
M = [("FC_cplex", "FasterCC", "CPLEX"), ("FASTCC_cplex", "FASTCC", "CPLEX"), ("SPECTRA_cplex", "SPECTRA", "CPLEX"),
     ("oldFC_cplex", "FasterCC (preprint)", "CPLEX"), ("FC_gurobi", "FasterCC", "Gurobi"), ("FASTCC_gurobi", "FASTCC", "Gurobi"), ("SPECTRA_gurobi", "SPECTRA", "Gurobi")]
rows = []
for k, m, s in M:
    have = [r for r in R if r[k] not in ("NaN", "")]
    rows.append([m, s, len({r["task"] for r in have}), sum(r["verdict"] == "BLOCKED" and r[k] == "1" for r in have),
                 sum(r["verdict"] == "CONSISTENT" and r[k] == "0" for r in have), sum(r["verdict"] == "BLOCKED" and r[k] == "0" for r in have),
                 sum(r["verdict"] == "CONSISTENT" and r[k] == "1" for r in have)])
nb = sum(r["verdict"] == "BLOCKED" for r in R); nc = len(R) - nb
table("tableR23_biomass", ["method", "solver", "tasks with a result (of 1,140)", "blocked biomass kept (rows)", "working biomass removed (rows)",
                           "blocked biomass removed", "working biomass kept"], rows, f"""
**Table R23 | Biomass reactions kept or removed wrongly, per method.** Every biomass reaction of every cluster task ({len(R):,} rows:
{nc:,} that can carry flux, {nb} blocked), decided by exact LPs in the task's input network (biomass check of 10 Nov,
fill-in 2). 'FasterCC (preprint)' = the earlier release run on 527 tasks. No method removes a biomass reaction that can
carry flux. The one blocked biomass FasterCC keeps with CPLEX is task CC_c02_l050_r09 (biomassreaction_2, kept by every CPLEX
method; the KEPTCHK case). Source: data/biocheck.""")
fig, ax = plt.subplots(figsize=(6.4, 2.8)); x = range(len(rows))
ax.bar(x, [r[3] for r in rows], color=[COL[r[0]] for r in rows], hatch=None, edgecolor="white")
for i, r in enumerate(rows):
    ax.text(i, r[3] + 0.3, str(r[3]), ha="center", fontsize=7)
ax.set_xticks(list(x)); ax.set_xticklabels([f"{r[0]}".replace(" (", "\n(") + f"\n{r[1]}" for r in rows], fontsize=6.5); ax.set_ylabel("blocked biomass reactions kept (rows)")
save(fig, "figR23_biomass_blocked_kept", f"""
**Figure R23 | Blocked biomass reactions left in the output, per method and solver.** Of {nb} biomass reactions that cannot carry
flux (exact LPs) in the 1,140 cluster tasks, the number each method keeps. No method removes a working biomass reaction
(0 in all rows of tableR23_biomass). Data: tableR23_biomass.csv; code: make_bio_wb_me.py.""")

# ---- SPECTRA on the whole-body models
L = [l.split() for l in open(os.path.join(src, "spectrawb_v1_table.txt")) if l.startswith(("Harvetta", "Harvey"))]
hdr = ["model", "eps", "mode", "rxns", "Crows", "inOK", "SPECTRA kept", "LPs", "failed LPs", "sec", "FVA consistent", "missed", "extra", "outOK", "rerun kept", "rerun LPs", "rerun failed"]
table("tableWB_spectra_wholebody", hdr, [[int(v) if v.isdigit() else v for v in l] for l in L], """
**Table WB | SPECTRA on the corrected whole-body models against FVA.** missed = FVA-consistent but removed by SPECTRA; extra =
kept by SPECTRA but FVA-blocked (|v| >= 0.99 eps on the FVA ranges); inOK/outOK = a steady state exists in the input / in
SPECTRA's output; rerun = SPECTRA on its own output (a stable method returns the same set). 'published' = SPECTRA as
published; 'withC' = with the coupling rows. Source: spectrawb_v1_table.txt (bench2026/spectrawb_out).""")
fig, ax = plt.subplots(figsize=(5.2, 2.6)); lab = [f"{l[0].split('_')[0]}\neps {l[1]}\n{l[2]}" for l in L]
ax.bar([i - 0.2 for i in range(len(L))], [int(l[12]) for l in L], 0.4, color=COL["SPECTRA"], label="kept although FVA-blocked (extra)")
ax.bar([i + 0.2 for i in range(len(L))], [int(l[11]) for l in L], 0.4, color="#999", label="removed although FVA-consistent (missed)")
ax.set_xticks(range(len(L))); ax.set_xticklabels(lab, fontsize=6); ax.set_ylabel("reactions"); ax.legend(frameon=False, fontsize=6.5)
save(fig, "figWB_spectra_wholebody", """
**Figure WB | SPECTRA's errors on the whole-body models (aafix) by tolerance and mode.** At eps 1e-4 SPECTRA equals FVA; at
eps 1e-6 it keeps 62-70 FVA-blocked reactions as published and 430-558 with the coupling rows. Data: tableWB_spectra_wholebody.csv.""")

# ---- ME model snapshot
T = list(csv.DictReader(open(os.path.join(src, "me_fva_vs_quad.tsv")), delimiter="\t"))
c = collections.Counter((r["fastercc_snap4"], r["fva_call"], r["spectra"]) for r in T)
name = {"A": "in A (proven consistent)", "J": "still to decide (J)", "neither": "removed earlier"}
rows = [[name[a], b, s, n] for (a, b, s), n in sorted(c.items(), key=lambda kv: -kv[1])]
table("tableME_snapshot", ["FasterCC after the singleton stage (snapshot 4)", "FVA call", "SPECTRA", "reactions"], rows, f"""
**Table ME | ME model (ME_dbl_E0_L1_00e_2, {len(T):,} reactions): FasterCC's state after its singleton stage against FVA and
SPECTRA.** FasterCC's run on this model had not finished (restart running); the table is its snapshot 4. FVA is a
floating-point FVA; the exact checks of 14-15 Nov (quad precision, SoPlex, Farkas) are in the report. Source: me_fva_vs_quad.tsv.""")
