"""check_biomass.py -- independent re-count of the biomass check (Table R23) from data/biocheck, plus the breakdowns
that show WHY each blocked biomass reaction is blocked.

Input : data/biocheck/bio_cNN_fill2.csv (one row per biomass reaction per cluster task, 1,140 tasks).
        verdict (CONSISTENT / BLOCKED) and maxflux come from exact LPs on the task's input network (biomass check of
        10 Nov, fill-in 2). Method columns: 1 = the method kept the reaction, 0 = removed it, NaN = no result.
Output: printed to the screen and written to out/check_biomass_output.txt (refuses to overwrite).
        Stops with an error if the re-count differs from out/tableR23_biomass.csv.

Run: python3 check_biomass.py   (standard library only)
"""
import collections, csv, glob, os, sys

here = os.path.dirname(os.path.abspath(__file__))
R = []
for f in sorted(glob.glob(os.path.join(here, "data", "biocheck", "bio_c*_fill2.csv"))):
    R += list(csv.DictReader(open(f)))
if not R:
    sys.exit("no data/biocheck/bio_c*_fill2.csv found")
out_path = os.path.join(here, "out", "check_biomass_output.txt")
if os.path.exists(out_path):
    sys.exit(f"{out_path} exists; remove it first")
lines = []
def p(s=""):
    print(s); lines.append(s)

blocked = [r for r in R if r["verdict"] == "BLOCKED"]
working = [r for r in R if r["verdict"] == "CONSISTENT"]
other = [r for r in R if r["verdict"] not in ("BLOCKED", "CONSISTENT") or r["status"] != "ok"]
p(f"rows: {len(R):,}  tasks: {len({r['task'] for r in R}):,}  working (CONSISTENT): {len(working):,}  blocked: {len(blocked)}"
  f"  other verdict or status not ok: {len(other)}")
p("biomass reactions per type: " + ", ".join(f"{k} {v:,}" for k, v in
  sorted(collections.Counter(r['reaction'].rsplit('_', 1)[0] for r in R).items())))
mf = [float(r["maxflux"]) for r in working]
p(f"largest flux of the working biomass reactions (exact LP): min {min(mf):g}, max {max(mf):g}")

p("\n1. Why the blocked biomass reactions are blocked")
closed = [r for r in blocked if float(r["lb"]) == 0 and float(r["ub"]) == 0]
openb = [r for r in blocked if r not in closed]
p(f"   closed by its own bounds in the task (lb = ub = 0, one of the reactions the task removed): {len(closed)}")
p(f"   bounds open (0..1000) but blocked by the reactions removed around it:                    {len(openb)}")
p("   by number of reactions removed in the task (level): " + ", ".join(
  f"{k}: {v}" for k, v in sorted(collections.Counter(r['level'] for r in blocked).items())))
p("   (no blocked biomass reaction at levels 000 and 010)" if not any(r['level'] in ('000', '010') for r in blocked) else "")

p("\n2. Table R23 re-counted (rows; tasks with a result = tasks where the method column is 0 or 1)")
cols = [("FasterCC", "CPLEX", "FC_cplex"), ("FASTCC", "CPLEX", "FASTCC_cplex"), ("SPECTRA", "CPLEX", "SPECTRA_cplex"),
        ("FasterCC (preprint)", "CPLEX", "oldFC_cplex"), ("FasterCC", "Gurobi", "FC_gurobi"),
        ("FASTCC", "Gurobi", "FASTCC_gurobi"), ("SPECTRA", "Gurobi", "SPECTRA_gurobi")]
p(f"   {'method':<20} {'solver':<7} {'tasks':>6} {'blocked kept':>13} {'working removed':>16} {'blocked removed':>16} {'working kept':>13}")
mine = []
for name, solver, c in cols:
    has = [r for r in R if r[c] in ("0", "1")]
    t = len({r["task"] for r in has})
    bk = sum(r[c] == "1" for r in has if r["verdict"] == "BLOCKED")
    wr = sum(r[c] == "0" for r in has if r["verdict"] == "CONSISTENT")
    br = sum(r[c] == "0" for r in has if r["verdict"] == "BLOCKED")
    wk = sum(r[c] == "1" for r in has if r["verdict"] == "CONSISTENT")
    mine.append([name, solver, t, bk, wr, br, wk])
    p(f"   {name:<20} {solver:<7} {t:>6,} {bk:>13} {wr:>16} {br:>16} {wk:>13,}")

p("\n3. Every blocked biomass reaction a method kept (task, reaction, bounds, why blocked)")
for name, solver, c in cols:
    for r in blocked:
        if r[c] == "1":
            why = "closed by bounds" if r in closed else "bounds open, blocked by the network"
            p(f"   {name} {solver}: {r['task']} {r['reaction']} [{r['lb']}, {r['ub']}] {why}")

ref = os.path.join(here, "out", "tableR23_biomass.csv")
if os.path.exists(ref):
    T = list(csv.reader(open(ref)))[1:]
    ok = len(T) == len(mine) and all(
        [a[0], a[1]] + [int(x.replace(",", "")) for x in a[2:]] == m for a, m in zip(T, mine))
    p(f"\nre-count equals out/tableR23_biomass.csv: {'yes' if ok else 'NO'}")
    if not ok:
        open(out_path, "w").write("\n".join(lines) + "\n"); sys.exit("re-count differs from tableR23_biomass.csv")
open(out_path, "w").write("\n".join(lines) + "\n")
print(f"\nwritten: {out_path}")
