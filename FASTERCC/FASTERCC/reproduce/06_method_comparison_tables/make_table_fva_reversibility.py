#!/usr/bin/env python3
"""make_table_fva_reversibility.py -- Table 3 extended to every model: reversible reactions by FVA against FasterCC.
For each model: FVA classes of the reactions (reversible, forward only, backward only, blocked; threshold 1e-6 as in Table 3),
reversible reactions by the input bounds, reversible reactions FasterCC's preprocessing made irreversible and whether FVA
agrees with the direction (reaction by reaction), reversible reactions left in FasterCC's consistent network, and the gap to
FVA: the reactions FasterCC leaves reversible although FVA says they can run in one direction only. The gap is what
FasterCC's optional reversibility check (5th argument revopt = 1) is meant to close; that check was NOT run with the final
version, so the gap is reported, not corrected.
Inputs (this folder): fcc_steps_EXA2_16Nov_models.csv and fcc_steps_EXA2_16Nov_reoriented.csv (FasterCC step records, default
mode, the runs used for Table 3), methods_real_16Nov_models.csv (final benchmark), tasks_real.tsv (FVA tag of each model),
fva_lbub/<tag>_FVA_lbub.csv (FVA minimum and maximum per reaction, columns rxn_id, lb, ub) and fva_lbub/FVA_iMM*.xlsx.
A model without an FVA file is listed with its FasterCC columns and 'FVA file not available'.
Usage: python3 -I make_table_fva_reversibility.py OUT.csv"""
import csv, sys, os, glob
if len(sys.argv) != 2: raise SystemExit(__doc__)
if os.path.exists(sys.argv[1]): raise SystemExit(sys.argv[1] + " exists: not overwriting")
EPS = 1e-6
ST = {m["model"]: m for m in csv.DictReader(open("fcc_steps_EXA2_16Nov_models.csv"))}
RO = {}
for r in csv.DictReader(open("fcc_steps_EXA2_16Nov_reoriented.csv")):
    if r["mode"] == "def": RO.setdefault(r["model"], []).append(r)
FN = {m["model"]: m for m in csv.DictReader(open("methods_real_16Nov_models.csv"))}
TAG = {}
for t in csv.DictReader(open("tasks_real.tsv"), delimiter="\t"): TAG.setdefault(t["model"], t["fva"])
G = {r["model"]: r["group"] for r in csv.DictReader(open("../02_FasterCC_vs_SPECTRA/out/figSP3_real_model_times.csv"))}
def fva(model):
    for tag in (model, TAG.get(model, "")):
        f = os.path.join("fva_lbub", f"{tag}_FVA_lbub.csv")
        if tag and os.path.exists(f):
            return {r["rxn_id"]: (float(r["lb"]), float(r["ub"])) for r in csv.DictReader(open(f))}, os.path.basename(f)
    return None, ""
rows = []
order = ["Small test", "AGORA microbes", "Curated microbial", "Human generic", "Fetal tissues", "Whole-body", "Microbiota", "ME model"]
models = sorted(FN, key=lambda m: (order.index(G[m]), int(FN[m]["in_total"]))) + sorted(m for m in ST if m.startswith("Cluster"))
for m in models:
    F, src = fva(m); s = ST.get(m); f = FN.get(m)
    r = dict(model=m, type=G.get(m, "Cluster network (no reaction closed)"), reactions=(f or s)["in_total"],
             fva_consistent=f["NEW_nFVA"] if f else "", fva_file=src or "FVA file not available")
    if F:
        cls = {k: ("rev" if a < -EPS and b > EPS else "fwd" if b > EPS else "bwd" if a < -EPS else "blocked") for k, (a, b) in F.items()}
        for c in ("rev", "fwd", "bwd", "blocked"): r["fva_" + {"rev": "reversible", "fwd": "forward_only", "bwd": "backward_only", "blocked": "blocked"}[c]] = sum(v == c for v in cls.values())
    else:
        for c in ("reversible", "forward_only", "backward_only", "blocked"): r["fva_" + c] = ""
    r["reversible_by_bounds"] = (f or s)["in_rev"]
    if s and m in RO:
        toirr = [x for x in RO[m] if x["input_class"] == "rev" and x["lp_class"] != "rev"]
        r["made_irreversible_by_preprocessing"] = len(toirr)
        r["reversible_in_fasterCC_network"] = s["outLP_rev"]
        if F:
            miss = [x for x in toirr if x["rxn"] not in F]
            same = sum((F[x["rxn"]][1] <= EPS) if x["flipped"] == "1" else (F[x["rxn"]][0] >= -EPS) for x in toirr if x["rxn"] in F)
            r["same_direction_in_fva"] = same; r["other_direction_in_fva"] = len(toirr) - len(miss) - same; r["ids_not_in_fva"] = len(miss)
            r["gap_left_reversible_but_one_direction_in_fva"] = int(s["outLP_rev"]) - r["fva_reversible"]
        else:
            r["same_direction_in_fva"] = r["other_direction_in_fva"] = r["ids_not_in_fva"] = r["gap_left_reversible_but_one_direction_in_fva"] = ""
    else:
        r["made_irreversible_by_preprocessing"] = r["reversible_in_fasterCC_network"] = "no step record"
        r["same_direction_in_fva"] = r["other_direction_in_fva"] = r["ids_not_in_fva"] = r["gap_left_reversible_but_one_direction_in_fva"] = ""
    r["reversible_in_final_benchmark_output"] = (f["NEWout_rev"] if f and f["NEW_finished"] == "1" else ("run not finished" if f else ""))
    r["optional_reversibility_check"] = "not run with the final version"
    rows.append(r)
import openpyxl
for name, fn in (("iMM1415", "FVA_iMM1415_fullmodel.xlsx"), ("iMM1865", "FVA_iMM1865_fullmodel.xlsx")):
    p = os.path.join("fva_lbub", fn)
    if not os.path.exists(p): continue
    R = list(openpyxl.load_workbook(p, read_only=True)["full_model"].iter_rows(min_row=2, values_only=True))
    r = {k: "" for k in rows[0]}; r.update(model=name, type="Preprint Table 3 (not in the benchmark)", reactions=len(R), fva_file=fn,
        reversible_by_bounds=sum(1 for x in R if x[4] < 0 < x[5]), made_irreversible_by_preprocessing="not run with the final version",
        reversible_in_fasterCC_network="not run with the final version", optional_reversibility_check="not run with the final version")
    for c, cond in (("reversible", lambda a, b: a < -EPS and b > EPS), ("forward_only", lambda a, b: not a < -EPS and b > EPS),
                    ("backward_only", lambda a, b: a < -EPS and not b > EPS), ("blocked", lambda a, b: not a < -EPS and not b > EPS)):
        r["fva_" + c] = sum(1 for x in R if cond(x[7], x[8]))
    rows.append(r)
cols = list(rows[0])
w = csv.DictWriter(open(sys.argv[1], "w", newline=""), fieldnames=cols); w.writeheader(); w.writerows(rows)
print(len(rows), "rows;", sum(1 for r in rows if r["fva_file"] != "FVA file not available"), "with an FVA file")
