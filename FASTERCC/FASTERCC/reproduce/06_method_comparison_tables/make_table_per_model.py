#!/usr/bin/env python3
"""make_table_per_model.py -- readable per-model table from methods_real_16Nov_models.csv (+ model type from
../02_FasterCC_vs_SPECTRA/out/figSP3_real_model_times.csv). Usage: python3 -I make_table_per_model.py OUT.csv"""
import csv, sys, math
M = list(csv.DictReader(open("methods_real_16Nov_models.csv")))
G = {r["model"]: r["group"] for r in csv.DictReader(open("../02_FasterCC_vs_SPECTRA/out/figSP3_real_model_times.csv"))}
def num(x):
    try: v = float(x); return None if math.isnan(v) else v
    except Exception: return None
fi = lambda x: "" if num(x) is None else f"{int(round(num(x))):,}"
ft = lambda x: "no result" if num(x) is None else (f"{num(x):,.0f}" if num(x) >= 100 else f"{num(x):.1f}")
pct = lambda a, b: "" if num(a) is None or not num(b) else f"{100 * num(a) / num(b):.1f}"
NOTE = {"Microbiota": "every method’s output has 1 reaction ID that is not in the input model",
        "ME_dbl_E0_L1_00e_2": "FasterCC run not finished; FASTCC and SPECTRA outputs have 6 reaction IDs that are not in the input model"}
order = ["Small test", "AGORA microbes", "Curated microbial", "Human generic", "Fetal tissues", "Whole-body", "Microbiota", "ME model"]
rows = []
for m in M:
    fin = {k: m[f"{k}_finished"] == "1" for k in ("NEW", "FASTCC", "SPECTRA")}
    v = lambda k, c, f=fi: f(m[f"{k}_{c}"]) if fin[k] else "no result"
    rows.append(dict(model=m["model"], type=G.get(m["model"], ""), reactions=fi(m["in_total"]), fva_consistent=fi(m["NEW_nFVA"]),
        kept_FasterCC=v("NEW", "n"), kept_FASTCC=v("FASTCC", "n"), kept_SPECTRA=v("SPECTRA", "n"),
        removed_consistent_FasterCC=v("NEW", "falseBlocked"), removed_consistent_FASTCC=v("FASTCC", "falseBlocked"), removed_consistent_SPECTRA=v("SPECTRA", "falseBlocked"),
        kept_blocked_FasterCC=v("NEW", "falseConsistent"), kept_blocked_FASTCC=v("FASTCC", "falseConsistent"), kept_blocked_SPECTRA=v("SPECTRA", "falseConsistent"),
        time_s_FasterCC=v("NEW", "time_s", ft), time_s_FASTCC=v("FASTCC", "time_s", ft), time_s_SPECTRA=v("SPECTRA", "time_s", ft),
        reversible_pct_FasterCC=pct(m["NEWout_rev"], m["NEW_n"]) if fin["NEW"] else "no result",
        reversible_pct_FASTCC=pct(m["FASTCC_rev"], m["FASTCC_n"]), reversible_pct_SPECTRA=pct(m["SPECTRA_rev"], m["SPECTRA_n"]),
        notes=NOTE.get(m["model"], "") if m["checks"] != "OK" else ""))
rows.sort(key=lambda r: (order.index(r["type"]) if r["type"] in order else 99, num(r["reactions"].replace(",", "")) or 0))
w = csv.DictWriter(open(sys.argv[1], "w", newline=""), fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(len(rows), "models written")
