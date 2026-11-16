#!/usr/bin/env python3
"""make_table_models_preprocessing.py -- per-model description and FasterCC preprocessing, step by step.
Reads (all in this folder): methods_real_16Nov_models.csv (final benchmark: input counts, kept, FVA), fcc_steps_EXA2_16Nov_models.csv
and fcc_steps_EXA2_16Nov_runs.csv (FasterCC step records: dead ends, Poolman, LP stage, reorientation; written by
fcc_step_table.py from the runs' own REPORT files), real_runs_16Nov.tsv (input conversion done before the methods),
tasks_real.tsv (epsilon of each model in the final benchmark), ../02_FasterCC_vs_SPECTRA/out/figSP3_real_model_times.csv (model type).
The step records come from FasterCC runs of build EXA2 (jobs 16199964-16200159, [FIX] version C2OBJ_FIX, auto mode), not from the
final benchmark runs; column kept_same_as_final checks that they keep the same number of reactions as the final benchmark.
Usage: python3 -I make_table_models_preprocessing.py OUT.csv"""
import csv, sys, os
if len(sys.argv) != 2: raise SystemExit(__doc__)
if os.path.exists(sys.argv[1]): raise SystemExit(sys.argv[1] + " exists: not overwriting")
F = {m["model"]: m for m in csv.DictReader(open("methods_real_16Nov_models.csv"))}
ST = {m["model"]: m for m in csv.DictReader(open("fcc_steps_EXA2_16Nov_models.csv"))}
RU = {}
for r in csv.DictReader(open("fcc_steps_EXA2_16Nov_runs.csv")):
    if r["mode"] == "auto": RU.setdefault(r["model"], r)
G = {r["model"]: r["group"] for r in csv.DictReader(open("../02_FasterCC_vs_SPECTRA/out/figSP3_real_model_times.csv"))}
EPS = {}
for t in csv.DictReader(open("tasks_real.tsv"), delimiter="\t"):
    if t["method"] == "NEW": EPS[t["model"]] = t["eps"]
NORM = {}
for r in csv.DictReader(open("real_runs_16Nov.tsv"), delimiter="\t"):
    if r["method"] == "NEW": NORM[r["model"]] = r["normalization"]
def plain(n):
    if n.startswith("already in FasterCC form"): return "none (already in the form S·v = 0)"
    if n == "csense missing -> E": return "missing constraint senses set to equality"
    if n == "b missing -> 0 | csense missing -> E | c missing -> 0": return "missing b, csense and c set to 0, equality and 0"
    if n.startswith("A/b/csense -> S (1668 rows)"): return "A matrix written as S (1,668 rows); no coupling rows"
    if n.startswith("duplicate ID EX_biomass[c]"): return "1 duplicate reaction ID renamed; A matrix split into S (154,801 rows) and coupling constraints C (245,811 rows)"
    if n.startswith("duplicate ID DM_leuW_tRNA"): return "6 duplicate reaction IDs renamed; 6,086 rows with a non-zero right-hand side or an inequality moved to coupling constraints C"
    return n
def eps(x):
    return {"0.0001": "1e-4", "1e-4": "1e-4", "1e-6": "1e-6", "1": "1"}.get(x, x)
order = ["Small test", "AGORA microbes", "Curated microbial", "Human generic", "Fetal tissues", "Whole-body", "Microbiota", "ME model"]
NR = "not recorded"
rows = []
def add(model, typ):
    f = F.get(model); s = ST.get(model); r = RU.get(model)
    lp = lambda: str(int(s["LP_removed_auto"] or 0))
    row = dict(model=model, type=typ,
        reactions=f["in_total"] if f else s["in_total"],
        metabolites=r["mets"] if r else NR,
        epsilon=eps(EPS.get(model, r["epsilon"] if r else "")),
        input_reversible=f["in_rev"] if f else s["in_rev"],
        input_irreversible_forward=f["in_fwd"] if f else s["in_irrev_fwd"],
        input_irreversible_backward=f["in_bwd"] if f else s["in_irrev_bwd"],
        input_closed=f["in_closed"] if f else s["in_closed"],
        removed_dead_ends=s["P1P2_deadend_initial"] if s else NR,
        removed_poolman_blocked=s["P3_poolman_blocked"] if s else NR,
        removed_dead_ends_after_poolman=s["P3_deadend_after_poolman"] if s else NR,
        removed_lp_stage=lp() if s else NR,
        kept_consistent=s["consistent_auto"] if s else (f["NEW_n"] if f and f["NEW_finished"] == "1" else "run not finished"),
        reversible_made_irreversible=str(int(s["rev_to_fwd_flipped"]) + int(s["rev_to_fwd_bounds"])) if s else NR,
        of_which_backward_only_flipped=s["rev_to_fwd_flipped"] if s else NR,
        of_which_forward_only_bounds=s["rev_to_fwd_bounds"] if s else NR,
        backward_irreversible_written_forward=s["fixIrr_bwd_to_fwd"] if s else NR,
        reversible_flipped_still_reversible=s["rev_flipped_still_rev"] if s else NR,
        columns_flipped_total=s["reoriented_at_LP"] if s else NR,
        reversible_at_lp=s["lp_rev"] if s else NR,
        reversible_in_output_step_record=s["outLP_rev"] if s else NR,
        reversible_in_output_final=(f["NEWout_rev"] if f and f["NEW_finished"] == "1" else ("run not finished" if f else "")),
        kept_final_benchmark=(f["NEW_n"] if f and f["NEW_finished"] == "1" else ("run not finished" if f else "")),
        kept_same_as_final=("yes" if s and f and f["NEW_finished"] == "1" and s["consistent_auto"] == f["NEW_n"] else ("" if not (s and f) else "no")),
        fva_consistent=f["NEW_nFVA"] if f else "",
        fasterCC_kept_fva_blocked_final=(f["NEW_falseConsistent"] if f and f["NEW_finished"] == "1" else ""),
        input_conversion=plain(NORM.get(model, "")) if f else "",
        step_record=("FasterCC build EXA2, auto mode, job " + r["job"]) if r else "none in the files available")
    if s:
        tot = sum(int(s[k] or 0) for k in ["P1P2_deadend_initial", "P3_poolman_blocked", "P3_deadend_after_poolman"]) + int(lp()) + int(s["consistent_auto"])
        if str(tot) != row["reactions"]: raise SystemExit("steps do not add up for " + model)
    rows.append(row)
for m in sorted(F, key=lambda m: (order.index(G[m]), int(F[m]["in_total"]))): add(m, G[m])
for m in sorted([m for m in ST if m.startswith("Cluster")]): add(m, "Cluster network (no reaction closed)")
w = csv.DictWriter(open(sys.argv[1], "w", newline=""), fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(len(rows), "rows written")
