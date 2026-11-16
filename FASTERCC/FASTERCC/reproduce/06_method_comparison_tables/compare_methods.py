#!/usr/bin/env python3
"""compare_methods.py -- READ ONLY. Comparison of FasterCC (NEW), FASTCC and SPECTRA from the TSV files written by
export_compare.m. No recomputation: every number is counted from those files.
Direction from the bounds: rev (lb<0<ub), fwd (lb>=0, ub>0), bwd (lb<0, ub<=0), closed (lb=ub=0); irrev = fwd + bwd.

Real models (real_rxns.tsv + real_runs.tsv) -> PREFIX_real_models.csv, one row per model:
  input total / rev / irrev (fwd, bwd) / closed;
  per method: finished, kept, kept rev / irrev by the INPUT bounds (what the saved model of each method contains);
  FasterCC also by its returned bounds (resModel: with its reorientation);
  kept by all three, only-FasterCC, only-FASTCC, only-SPECTRA, FasterCC vs FASTCC / SPECTRA differences;
  FVA check of each method (nFVA, falseBlocked, falseConsistent) as stored by the benchmark.
  PREFIX_real_diff_ids.csv: every reaction kept by some but not all finished methods, with its bounds.
Benchmark (bench_tasks.tsv + bench_diffs.tsv) -> PREFIX_bench_tasks.csv (as exported, plus irrev counts) and
  PREFIX_bench_summary.csv: per folder, cluster and method: tasks, kept and in rev/irrev (medians), tasks where the
  three sets are identical, and how many tasks had any difference; PREFIX_bench_diff_ids.csv: the differing IDs.
Checks printed: input rev count = res.nRev in every task; kept counts per method = n; models with no result listed.
Usage: python3 -I compare_methods.py DIR --out PREFIX"""
import csv, os, sys, statistics, collections
a = sys.argv[1:]; out = a[a.index("--out") + 1]; D = [x for x in a if x not in ("--out", out)][0]
for s in ("_real_models.csv", "_real_diff_ids.csv", "_bench_tasks.csv", "_bench_summary.csv", "_bench_diff_ids.csv"):
    if os.path.exists(out + s): raise SystemExit(out + s + " exists: not overwriting")
def cls(lb, ub):
    try: lb, ub = float(lb), float(ub)
    except ValueError: return "NA"
    if lb == 0 and ub == 0: return "closed"
    return "rev" if lb < 0 < ub else "fwd" if (lb >= 0 and ub > 0) else "bwd" if (lb < 0 and ub <= 0) else "other"
M = ["NEW", "FASTCC", "SPECTRA"]
# ---------------- real models
runs = collections.defaultdict(dict)
for r in csv.DictReader(open(os.path.join(D, "real_runs.tsv")), delimiter="\t"): runs[r["model"]][r["method"]] = r
C = collections.defaultdict(collections.Counter); diff_rows = []
for r in csv.DictReader(open(os.path.join(D, "real_rxns.tsv")), delimiter="\t"):
    m = r["model"]; c = C[m]; ci = cls(r["lb_in"], r["ub_in"]); c["in_total"] += 1; c["in_" + ci] += 1
    k = {x: r["kept_" + x] for x in M}; fin = [x for x in M if k[x] != "NA"]
    for x in fin:
        if k[x] == "1": c[x + "_n"] += 1; c[x + "_in_" + ci] += 1
    if k["NEW"] == "1": c["NEWout_" + cls(r["lb_NEWout"], r["ub_NEWout"])] += 1
    kk = {x: k[x] == "1" for x in fin}
    if len(fin) == 3 and all(kk.values()): c["all3"] += 1
    if len(fin) >= 2 and len(set(kk.values())) > 1:
        tag = "only_" + "+".join(x for x in fin if kk[x])
        c[tag] += 1
        diff_rows.append(dict(model=m, rxn=r["rxn"], kept_by="+".join(x for x in fin if kk[x]),
                              not_kept_by="+".join(x for x in fin if not kk[x]), lb_in=r["lb_in"], ub_in=r["ub_in"],
                              input_class=ci, lb_NEWout=r["lb_NEWout"], ub_NEWout=r["ub_NEWout"]))
cols = ["model", "in_total", "in_rev", "in_irrev", "in_fwd", "in_bwd", "in_closed"]
for x in M: cols += [f"{x}_finished", f"{x}_n", f"{x}_rev", f"{x}_irrev", f"{x}_nFVA", f"{x}_falseBlocked", f"{x}_falseConsistent", f"{x}_time_s"]
cols += ["NEWout_rev", "NEWout_irrev", "all3", "only_NEW", "only_FASTCC", "only_SPECTRA", "NEW+FASTCC_not_SPECTRA",
         "NEW+SPECTRA_not_FASTCC", "FASTCC+SPECTRA_not_NEW", "checks"]
rows = []
for m in sorted(C):
    c = C[m]; R = runs.get(m, {}); row = dict(model=m, in_total=c["in_total"], in_rev=c["in_rev"], in_irrev=c["in_fwd"] + c["in_bwd"],
                                               in_fwd=c["in_fwd"], in_bwd=c["in_bwd"], in_closed=c["in_closed"])
    bad = []
    for x in M:
        r = R.get(x, {}); fin = r.get("finished", "NA")
        row.update({f"{x}_finished": fin, f"{x}_n": c[x + "_n"] if fin == "1" else "", f"{x}_rev": c[x + "_in_rev"] if fin == "1" else "",
                    f"{x}_irrev": c[x + "_in_fwd"] + c[x + "_in_bwd"] if fin == "1" else "", f"{x}_nFVA": r.get("nFVA", ""),
                    f"{x}_falseBlocked": r.get("nFalseBlocked", ""), f"{x}_falseConsistent": r.get("nFalseConsistent", ""),
                    f"{x}_time_s": r.get("time_s", "")})
        if fin == "1" and str(c[x + "_n"]) != r.get("n"): bad.append(f"{x} kept {c[x + '_n']} != n {r.get('n')}")
        if fin == "1" and r.get("ids_not_in_input", "0") not in ("0", ""): bad.append(f"{x} {r['ids_not_in_input']} ids not in input")
    row.update(NEWout_rev=c["NEWout_rev"], NEWout_irrev=c["NEWout_fwd"] + c["NEWout_bwd"], all3=c["all3"],
               only_NEW=c["only_NEW"], only_FASTCC=c["only_FASTCC"], only_SPECTRA=c["only_SPECTRA"],
               **{"NEW+FASTCC_not_SPECTRA": c["only_NEW+FASTCC"], "NEW+SPECTRA_not_FASTCC": c["only_NEW+SPECTRA"],
                  "FASTCC+SPECTRA_not_NEW": c["only_FASTCC+SPECTRA"]}, checks="OK" if not bad else " | ".join(bad))
    rows.append(row)
with open(out + "_real_models.csv", "x", newline="") as fh: w = csv.DictWriter(fh, cols); w.writeheader(); w.writerows(rows)
with open(out + "_real_diff_ids.csv", "x", newline="") as fh:
    w = csv.DictWriter(fh, list(diff_rows[0]) if diff_rows else ["model"]); w.writeheader(); w.writerows(diff_rows)
print(f"real: {len(rows)} models, checks failed: {sum(r['checks'] != 'OK' for r in rows)}, differing reaction rows: {len(diff_rows)}")
for x in M: print(f"  {x}: finished {sum(r[f'{x}_finished'] == '1' for r in rows)}, no result {[r['model'] for r in rows if r[f'{x}_finished'] != '1']}")
# ---------------- benchmark
T = list(csv.DictReader(open(os.path.join(D, "bench_tasks.tsv")), delimiter="\t"))
for t in T: t["kept_irrev"] = int(t["kept_fwd"]) + int(t["kept_bwd"]); t["in_irrev"] = int(t["in_fwd"]) + int(t["in_bwd"])
with open(out + "_bench_tasks.csv", "x", newline="") as fh: w = csv.DictWriter(fh, list(T[0])); w.writeheader(); w.writerows(T)
Dd = collections.defaultdict(list)
for d in csv.DictReader(open(os.path.join(D, "bench_diffs.tsv")), delimiter="\t"): Dd[(d["folder"], d["task"])].append(d)
with open(out + "_bench_diff_ids.csv", "x", newline="") as fh:
    w = csv.writer(fh); w.writerow(["folder", "task", "methodA", "methodB", "side", "rxn"])
    for k in sorted(Dd):
        for d in Dd[k]: w.writerow([d["folder"], d["task"], d["methodA"], d["methodB"], d["side"], d["rxn"]])
G = collections.defaultdict(list)
for t in T: G[(t["folder"], int(t["cluster"]), t["method"])].append(t)
nset = collections.defaultdict(set)
for t in T: nset[(t["folder"], t["task"])].add(t["n"])
with open(out + "_bench_summary.csv", "x", newline="") as fh:
    w = csv.writer(fh); w.writerow(["folder", "cluster", "method", "tasks", "median_in_total", "median_in_rev", "median_in_irrev",
                                    "median_kept", "median_kept_rev", "median_kept_irrev", "median_nRev_reported",
                                    "tasks_with_id_differences_in_file", "median_time_s"])
    for k in sorted(G):
        L = G[k]; med = lambda f: statistics.median(float(x[f]) for x in L)
        nr = [float(x["nRev_reported"]) for x in L if x["nRev_reported"] not in ("NA", "NaN", "")]
        w.writerow([k[0], k[1], k[2], len(L), med("in_total"), med("in_rev"), med("in_irrev"), med("n"), med("kept_rev"),
                    med("kept_irrev"), statistics.median(nr) if nr else "", sum(1 for x in L if Dd.get((x["folder"], x["task"]))),
                    med("time_s")])
bad = [t for t in T if t["in_rev_check"] != "OK"]
print(f"bench: {len(T)} rows, {len({(t['folder'], t['task']) for t in T})} task files; input rev check failed: {len(bad)}; "
      f"task files with any ID difference: {len(Dd)}; per folder/method rows: " +
      ", ".join(f"{f}/{m} {n}" for (f, m), n in sorted(collections.Counter((t['folder'], t['method']) for t in T).items())))
print("written", ", ".join(out + s for s in ("_real_models.csv", "_real_diff_ids.csv", "_bench_tasks.csv", "_bench_summary.csv", "_bench_diff_ids.csv")))
