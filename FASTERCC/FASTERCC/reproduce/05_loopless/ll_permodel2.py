#!/usr/bin/env python3
"""ll_permodel2.py -- READ ONLY. Copy of ll_permodel.py (md5 15c4bee4) plus --csv PREFIX: writes PREFIX_permodel.csv
(one row per model and method) and PREFIX_pertype.csv (one row per model type and method) for plots. Text output unchanged.
 Loopless results PER MODEL and PER MODEL TYPE, from loopkit/llcheck/ll_table.tsv and the
checker files loopkit/llcheck/<model>__<NET>_llcheck_all_s*of*.tsv (the same files ll_table.py reads).
Per model and method:
  capable      loop-capable reactions
  stepB        certified loop-free by step B (joint LP)
  island       loop-only by closed islands (S alone)
  checker      decided by the checker = (FINAL_loop_free - stepB) + (FINAL_loop_only - island)
  undec        FINAL_undecided, and undec % of capable
  tries        checker tries (candidate directions), timeTL = tries stopped at the time limit, and its %
Model type: assigned by name rules printed in the output (check them; anything not matched is 'other').
Usage: python3 -I ll_permodel2.py [--out FILE] [--csv PREFIX]   (refuses to overwrite any output file)"""
import csv, glob, os, re, sys, statistics, collections
LK = "/mnt/scratch/users/egonzalez/fastbox2copy/bench2026/loopkit"
out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
csvp = sys.argv[sys.argv.index("--csv") + 1] if "--csv" in sys.argv else None
fo = open(out, "x") if out else sys.stdout
def p(*a): print(*a, file=fo)
RULES = [("whole-body", r"^(Harvey|Harvetta)"), ("community", r"^Microbiota$"), ("fetal tissue", r"^Fetus_"),
         ("human generic", r"^(Recon|Human)"), ("mitochondrial", r"^cardiac_mit"),
         ("BiGG microbe", r"^(i[A-Z]|Ec_|Sc_|STM_)"), ("AGORA bacterium", r"^[A-Z][a-z]+_[a-z]+")]
def mtype(m):
    for t, rx in RULES:
        if re.search(rx, m): return t
    return "other"
R = list(csv.DictReader(open(f"{LK}/llcheck/ll_table.tsv"), delimiter="\t"))
TR = collections.defaultdict(collections.Counter)            # (model, method) -> verdict counts of candidate tries
for g in glob.glob(f"{LK}/llcheck/*_llcheck_all_s*of*.tsv"):
    net = os.path.basename(g).split("_llcheck_all_")[0]; model, _, meth = net.rpartition("__")
    meth = {"NEW": "FasterCC"}.get(meth, meth)
    for r in csv.DictReader(open(g), delimiter="\t"):
        if r["kind"] == "candidate": TR[(model, meth)][r["verdict"]] += 1
I = lambda r, c: int(r[c])
rows = []
for r in R:
    m, me = r["model"], r["method"]; c = I(r, "loop_capable")
    chk = (I(r, "FINAL_loop_free") - I(r, "B_loopfree_LP")) + (I(r, "FINAL_loop_only") - I(r, "island_looponly"))
    t = TR.get((m, me), collections.Counter()); nt = sum(t.values()); tl = t["UNDECIDED_time_limit"]
    rows.append(dict(model=m, type=mtype(m), method=me, capable=c, stepB=I(r, "B_loopfree_LP"), island=I(r, "island_looponly"),
                     checker=chk, undec=I(r, "FINAL_undecided"), undec_pct=(100 * I(r, "FINAL_undecided") / c if c else None),
                     tries=nt, timeTL=tl, tl_pct=(100 * tl / nt if nt else None), shards=r["checker_shards_done"]))
f1 = lambda x: "-" if x is None else f"{x:.1f}"
p("MODEL TYPE RULES (name -> type):", "; ".join(f"{t}: {rx}" for t, rx in RULES), "; anything else: other")
p("\n=== 1) PER MODEL (one line per method)")
p(f"{'type':16s} {'model':38s} {'method':9s} {'capable':>7s} {'stepB':>6s} {'island':>6s} {'checker':>7s} {'undec':>6s} {'undec%':>6s} {'tries':>6s} {'timeTL':>6s} {'TL%':>6s} shards")
for x in sorted(rows, key=lambda x: (x["type"], x["model"], x["method"])):
    p(f"{x['type']:16s} {x['model'][:38]:38s} {x['method']:9s} {x['capable']:7d} {x['stepB']:6d} {x['island']:6d} {x['checker']:7d} {x['undec']:6d} {f1(x['undec_pct']):>6s} {x['tries']:6d} {x['timeTL']:6d} {f1(x['tl_pct']):>6s} {x['shards']}")
p("\n=== 2) PER MODEL TYPE AND METHOD (per-model values; models with 0 loop-capable left out of the medians)")
p(f"{'type':16s} {'method':9s} {'models':>6s} {'stepB% median':>13s} {'undec% median':>13s} {'undec% min-max':>15s} {'TL% median':>10s} {'undec sum':>9s}")
G = collections.defaultdict(list)
for x in rows: G[(x["type"], x["method"])].append(x)
for (t, me), xs in sorted(G.items()):
    ys = [x for x in xs if x["capable"] > 0]
    sb = [100 * x["stepB"] / x["capable"] for x in ys]; up = [x["undec_pct"] for x in ys]; tp = [x["tl_pct"] for x in ys if x["tl_pct"] is not None]
    md = lambda v: f1(statistics.median(v)) if v else "-"
    rng = f"{min(up):.1f}-{max(up):.1f}" if up else "-"
    p(f"{t:16s} {me:9s} {len(xs):6d} {md(sb):>13s} {md(up):>13s} {rng:>15s} {md(tp):>10s} {sum(x['undec'] for x in xs):9d}")
p("\n=== 3) SAME MODEL, METHODS SIDE BY SIDE: undecided (undec%)  [INPUT = the model before any method]")
M = collections.defaultdict(dict)
for x in rows: M[x["model"]][x["method"]] = x
meths = ["INPUT", "FasterCC", "FASTCC", "SPECTRA"]
p(f"{'type':16s} {'model':38s} " + " ".join(f"{m:>16s}" for m in meths) + "   FasterCC-FASTCC FasterCC-SPECTRA")
for mdl in sorted(M, key=lambda k: (mtype(k), k)):
    d = M[mdl]; cell = lambda me: f"{d[me]['undec']} ({f1(d[me]['undec_pct'])})" if me in d else "no row"
    df = lambda a, b: f"{d[a]['undec'] - d[b]['undec']:+d}" if a in d and b in d else "-"
    p(f"{mtype(mdl):16s} {mdl[:38]:38s} " + " ".join(f"{cell(me):>16s}" for me in meths) + f"   {df('FasterCC','FASTCC'):>15s} {df('FasterCC','SPECTRA'):>16s}")
if csvp:
    K = ["type", "model", "method", "capable", "stepB", "island", "checker", "undec", "undec_pct", "tries", "timeTL", "tl_pct", "shards"]
    with open(csvp + "_permodel.csv", "x", newline="") as fh:
        w = csv.writer(fh); w.writerow(K)
        for x in sorted(rows, key=lambda x: (x["type"], x["model"], x["method"])):
            w.writerow(["" if x[k] is None else (round(x[k], 3) if isinstance(x[k], float) else x[k]) for k in K])
    with open(csvp + "_pertype.csv", "x", newline="") as fh:
        w = csv.writer(fh); w.writerow(["type", "method", "models", "models_with_loops", "stepB_pct_median", "undec_pct_median", "undec_pct_min", "undec_pct_max", "tl_pct_median", "undec_sum", "capable_sum"])
        for (t, me), xs in sorted(G.items()):
            ys = [x for x in xs if x["capable"] > 0]
            sb = [100 * x["stepB"] / x["capable"] for x in ys]; up = [x["undec_pct"] for x in ys]; tp = [x["tl_pct"] for x in ys if x["tl_pct"] is not None]
            r3 = lambda v: round(v, 3) if v is not None else ""
            w.writerow([t, me, len(xs), len(ys), r3(statistics.median(sb)) if sb else "", r3(statistics.median(up)) if up else "",
                        r3(min(up)) if up else "", r3(max(up)) if up else "", r3(statistics.median(tp)) if tp else "",
                        sum(x["undec"] for x in xs), sum(x["capable"] for x in xs)])
    print("written", csvp + "_permodel.csv", csvp + "_pertype.csv")
if out: fo.close(); print("written", out)
