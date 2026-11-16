#!/usr/bin/env python3
"""fcc_step_table.py -- READ ONLY. Per-model table of what FasterCC removed at each step, taken only from the files each
run wrote itself (no recomputation):
  REPORT/FATE_SUMMARY.csv      reactions per fate (dead ends P1/P2, Poolman P3, dead ends after Poolman, LP removals, consistent)
  REPORT/TIMELINE.csv          reactions left after each stage
  REPORT/REACTION_FATE.csv     per reaction: flippedAtLP (reoriented by the preprocessing), boundsChangedAtLP, forcedInInput
  CONSISTENT_RXNS.csv          the kept reactions (count must equal the CONSISTENT fates)
  RUN_LOG.txt                  model size, epsilon, version line, [AUTO] switches (def = FAST 0 COUPLED 0, auto = FAST 1
                               COUPLED 1), ALGORITHM TOTAL time, fccRestoreOutputModel "ALL CHECKS OK"
Optional --fva DIR: <model>_<def|auto>_diff.csv from fva_multi; rows with FVA=0, FCC=1 are the guard-kept reactions.
Checks per run (column 'checks'): fates sum to the input size; consistent fates = CONSISTENT_RXNS rows = final stage;
output-model checks OK. Any failed check is written, not hidden.
Usage: python3 -I fcc_step_table.py ROOT --out PREFIX [--fva DIR]   ROOT = folder holding <model>/RUN1_*/ run folders
Writes PREFIX_runs.csv (one row per run), PREFIX_models.csv (one row per model, def and auto side by side) and
PREFIX_reoriented.csv (every reaction whose direction or bounds were changed for the LP step, one run per model: def,
or auto if there is no def). Direction classes from the bounds: rev (lb<0<ub), fwd (lb>=0, ub>0), bwd (lb<0, ub<=0).
Counts per run (columns in_*, lp_*, out_*, outLP_*): total, reversible, irreversible (forward + backward), closed
(lb = ub = 0), for the input, the network at the LP, the output with the input bounds, and the output with the LP bounds.
Reorientation categories (input class -> class at the LP, flippedAtLP, boundsChangedAtLP):
  fixIrr_bwd_to_fwd        bwd -> fwd, flipped            backward-only irreversible written forward (fixIrr)
  rev_to_fwd_flipped       rev -> fwd, flipped, bounds    reversible that can only run backward: made forward-irreversible
  rev_to_fwd_bounds        rev -> fwd, not flipped, bounds reversible that can only run forward: made irreversible
  rev_flipped_still_rev    rev -> rev, flipped            reversible, sign of the column changed, still reversible
  other                    any other combination (listed with its classes)"""
import csv, glob, os, re, sys, collections
a = sys.argv[1:]
def opt(k):
    if k in a: v = a[a.index(k) + 1]; del a[a.index(k):a.index(k) + 2]; return v
out = opt("--out"); fva = opt("--fva"); root = a[0] if a else None
if not (out and root): raise SystemExit(__doc__)
for s in ("_runs.csv", "_models.csv", "_reoriented.csv"):
    if os.path.exists(out + s): raise SystemExit(out + s + " exists: not overwriting")
FATES = ["P1P2_deadend_initial", "P3_poolman_blocked", "P3_deadend_after_poolman", "LP_batch_removed",
         "LP_J_LP3_blocked", "LP_J_cert_skipped_blocked", "LP_A_removed_later", "CONSISTENT", "CONSISTENT_via_LP3"]
def typ(m):
    if m.startswith("Cluster"): return "cluster network"
    if m.startswith("Fetus_"): return "fetal tissue"
    if m.startswith(("Harvey", "Harvetta")): return "whole-body"
    if m.startswith(("Recon", "Human")): return "human generic"
    if m in ("cardiac_mit_glcuptake_atpmax",): return "mitochondrial"
    if m.startswith(("i", "Ec_", "Sc_", "STM", "ecoli")): return "BiGG microbe"
    if m in ("DAS", "subDas"): return "other"
    return "AGORA bacterium"
rows = []
for log in sorted(glob.glob(os.path.join(root, "*", "*", "RUN_LOG.txt"))):
    d = os.path.dirname(log); model = os.path.basename(os.path.dirname(d)); run = os.path.basename(d)
    r = dict(model=model, type=typ(model), run_folder=run, job=run.rsplit("_", 1)[-1])
    L = open(log, errors="replace").read()
    m = re.search(r"^\[AUTO\] FAST (\d) COUPLED (\d)", L, re.M)
    r["mode"] = {("0", "0"): "def", ("1", "1"): "auto"}.get(m.groups(), "FAST%s_COUPLED%s" % m.groups()) if m else "?"
    m = re.search(r"model: (\d+) reactions, (\d+) metabolites", L); r["input_rxns"], r["mets"] = (m.groups() if m else ("", ""))
    m = re.search(r"epsilon ([0-9.eE+-]+)\)", L); r["epsilon"] = m.group(1) if m else ""
    m = re.search(r"\[FIX\] version (\S+)", L); r["fix_version"] = m.group(1).rstrip(",") if m else ""
    m = re.search(r"ALGORITHM TOTAL.*?time\s+([0-9.]+) s", L); r["time_s"] = m.group(1) if m else ""
    r["output_checks"] = "OK" if re.search(r"ALL CHECKS\s+OK", L) else "NOT FOUND"
    F = collections.Counter()
    for x in csv.DictReader(open(os.path.join(d, "REPORT", "FATE_SUMMARY.csv"))): F[x["fate"]] += int(x["reactions"])
    for k in FATES: r[k] = F.get(k, 0)
    r["other_fates"] = ";".join(f"{k}={v}" for k, v in F.items() if k not in FATES)
    T = {x["stage"]: x["reactions"] for x in csv.DictReader(open(os.path.join(d, "REPORT", "TIMELINE.csv")))}
    r["after_P1P2"] = T.get("after P1 fixIrr + P2 dead ends", ""); r["after_P3"] = T.get("after P3 Poolman/dead-end closure", "")
    r["entering_LP7"] = T.get("after P4 dead ends = entering batch LP7", ""); r["final"] = T.get("final consistent", "")
    RF = list(csv.DictReader(open(os.path.join(d, "REPORT", "REACTION_FATE.csv"))))
    r["reoriented_at_LP"] = sum(x["flippedAtLP"] == "1" for x in RF)
    def dcls(lb, ub):
        try: lb, ub = float(lb), float(ub)
        except ValueError: return "NA"
        return "rev" if lb < 0 < ub else "fwd" if (lb >= 0 and ub > 0) else "bwd" if (lb < 0 and ub <= 0) else "zero"
    TR = []
    for x in RF:
        if x["presentAtLP"] != "1" or (x["flippedAtLP"] != "1" and x["boundsChangedAtLP"] != "1"): continue
        a, b, fl, bc = dcls(x["lbIn"], x["ubIn"]), dcls(x["lbAtLP"], x["ubAtLP"]), x["flippedAtLP"], x["boundsChangedAtLP"]
        cat = ("fixIrr_bwd_to_fwd" if (a, b, fl) == ("bwd", "fwd", "1") else
               "rev_to_fwd_flipped" if (a, b, fl, bc) == ("rev", "fwd", "1", "1") else
               "rev_to_fwd_bounds" if (a, b, fl, bc) == ("rev", "fwd", "0", "1") else
               "rev_flipped_still_rev" if (a, b, fl) == ("rev", "rev", "1") else f"other_{a}_to_{b}_flip{fl}_bounds{bc}")
        TR.append(dict(rxn=x["rxn"], input_class=a, lp_class=b, flipped=fl, bounds_changed=bc, category=cat, lbIn=x["lbIn"],
                       ubIn=x["ubIn"], lbAtLP=x["lbAtLP"], ubAtLP=x["ubAtLP"], fate=x["fate"], consistent=x["consistent"]))
    K = collections.Counter(t["category"] for t in TR)
    for k in ("fixIrr_bwd_to_fwd", "rev_to_fwd_flipped", "rev_to_fwd_bounds", "rev_flipped_still_rev"): r[k] = K.get(k, 0)
    r["reorient_other"] = sum(v for k, v in K.items() if k.startswith("other"))
    r["input_rev"] = sum(dcls(x["lbIn"], x["ubIn"]) == "rev" for x in RF); r["input_bwd"] = sum(dcls(x["lbIn"], x["ubIn"]) == "bwd" for x in RF)
    # reversible / irreversible counts at three points: input (all reactions), at the LP (present after preprocessing,
    # with the LP bounds), output (kept reactions; the output model keeps the input bounds; the LP-bound count shows
    # what the kept set would be with FasterCC's reorientation applied)
    def cnt(pre, sel, lbk, ubk):
        c = collections.Counter(dcls(x[lbk], x[ubk]) for x in RF if sel(x))
        r[pre + "_total"] = sum(c.values()); r[pre + "_rev"] = c["rev"]; r[pre + "_irrev"] = c["fwd"] + c["bwd"]
        r[pre + "_irrev_fwd"] = c["fwd"]; r[pre + "_irrev_bwd"] = c["bwd"]; r[pre + "_closed"] = c["zero"] + c["NA"]
    cnt("in", lambda x: True, "lbIn", "ubIn")
    cnt("lp", lambda x: x["presentAtLP"] == "1", "lbAtLP", "ubAtLP")
    cnt("out", lambda x: x["consistent"] == "1", "lbIn", "ubIn")
    cnt("outLP", lambda x: x["consistent"] == "1", "lbAtLP", "ubAtLP")
    r["_TR"] = TR
    r["bounds_changed_at_LP"] = sum(x["boundsChangedAtLP"] == "1" for x in RF)
    r["forced_in_input"] = sum(x["forcedInInput"] == "1" for x in RF)
    r["forced_not_consistent"] = sum(x["forcedInInput"] == "1" and x["consistent"] != "1" for x in RF)
    ncons = sum(1 for _ in csv.reader(open(os.path.join(d, "CONSISTENT_RXNS.csv")))) if os.path.exists(os.path.join(d, "CONSISTENT_RXNS.csv")) else -1
    r["consistent_total"] = F["CONSISTENT"] + F["CONSISTENT_via_LP3"]; r["removed_total"] = sum(F.values()) - r["consistent_total"]
    bad = []
    if str(sum(F.values())) != r["input_rxns"]: bad.append(f"fates {sum(F.values())} != input {r['input_rxns']}")
    if ncons not in (r["consistent_total"], r["consistent_total"] + 1): bad.append(f"CONSISTENT_RXNS rows {ncons} != {r['consistent_total']}")
    if r["final"] and int(r["final"]) != r["consistent_total"]: bad.append(f"final stage {r['final']} != {r['consistent_total']}")
    if r["output_checks"] != "OK": bad.append("output checks not found")
    r["checks"] = "OK" if not bad else " | ".join(bad)
    r["guard_kept_vs_FVA"] = ""
    if fva:
        f = os.path.join(fva, f"{model}_{r['mode']}_diff.csv")
        if os.path.exists(f):
            R = [x.rstrip("\n").split(",") for x in open(f)][1:]
            r["guard_kept_vs_FVA"] = sum(1 for x in R if len(x) >= 7 and x[-6] == "0" and x[-5] == "1")
    rows.append(r)
cols = ["model", "type", "mode", "job", "input_rxns", "mets", "epsilon", "fix_version", "P1P2_deadend_initial", "after_P1P2",
        "P3_poolman_blocked", "P3_deadend_after_poolman", "after_P3", "entering_LP7", "LP_batch_removed", "LP_J_LP3_blocked",
        "LP_J_cert_skipped_blocked", "LP_A_removed_later", "other_fates", "removed_total", "CONSISTENT", "CONSISTENT_via_LP3",
        "consistent_total", "final", "in_total", "in_rev", "in_irrev", "in_irrev_fwd", "in_irrev_bwd", "in_closed", "lp_total", "lp_rev", "lp_irrev", "lp_irrev_fwd", "lp_irrev_bwd", "lp_closed", "out_total", "out_rev", "out_irrev", "out_irrev_fwd", "out_irrev_bwd", "out_closed", "outLP_total", "outLP_rev", "outLP_irrev", "outLP_irrev_fwd", "outLP_irrev_bwd", "outLP_closed", "input_rev", "input_bwd", "reoriented_at_LP", "fixIrr_bwd_to_fwd",
        "rev_to_fwd_flipped", "rev_to_fwd_bounds", "rev_flipped_still_rev", "reorient_other", "bounds_changed_at_LP", "forced_in_input", "forced_not_consistent",
        "guard_kept_vs_FVA", "time_s", "output_checks", "checks", "run_folder"]
with open(out + "_runs.csv", "x", newline="") as fh:
    w = csv.DictWriter(fh, cols, extrasaction="ignore"); w.writeheader(); w.writerows(sorted(rows, key=lambda r: (r["type"], r["model"], r["mode"], r["job"])))
# per model: def and auto side by side; same-set flag
G = collections.defaultdict(dict)
for r in rows: G[r["model"]].setdefault(r["mode"], r)          # first run per mode (by job order) if a model has 4 runs
mc = ["model", "type", "input_rxns", "P1P2_deadend_initial", "P3_poolman_blocked", "P3_deadend_after_poolman", "LP_removed_def",
      "LP_removed_auto", "consistent_def", "consistent_auto", "in_total", "in_rev", "in_irrev", "in_irrev_fwd", "in_irrev_bwd", "in_closed", "lp_total", "lp_rev", "lp_irrev", "lp_irrev_fwd", "lp_irrev_bwd", "lp_closed", "out_total", "out_rev", "out_irrev", "out_irrev_fwd", "out_irrev_bwd", "out_closed", "outLP_total", "outLP_rev", "outLP_irrev", "outLP_irrev_fwd", "outLP_irrev_bwd", "outLP_closed", "input_rev", "input_bwd", "fixIrr_bwd_to_fwd",
      "rev_to_fwd_flipped", "rev_to_fwd_bounds", "rev_flipped_still_rev", "reorient_other", "reoriented_at_LP", "guard_kept_vs_FVA", "time_def_s", "time_auto_s", "runs", "checks"]
with open(out + "_models.csv", "x", newline="") as fh:
    w = csv.writer(fh); w.writerow(mc)
    for m in sorted(G, key=lambda m: (typ(m), m)):
        D, A = G[m].get("def", {}), G[m].get("auto", {}); B = D or A
        lp = lambda x: (x["LP_batch_removed"] + x["LP_J_LP3_blocked"] + x["LP_J_cert_skipped_blocked"] + x["LP_A_removed_later"]) if x else ""
        n = sum(1 for r in rows if r["model"] == m)
        w.writerow([m, B["type"], B["input_rxns"], B["P1P2_deadend_initial"], B["P3_poolman_blocked"], B["P3_deadend_after_poolman"],
                    lp(D), lp(A), D.get("consistent_total", ""), A.get("consistent_total", ""), *[B[k] for k in ["in_total", "in_rev", "in_irrev", "in_irrev_fwd", "in_irrev_bwd", "in_closed", "lp_total", "lp_rev", "lp_irrev", "lp_irrev_fwd", "lp_irrev_bwd", "lp_closed", "out_total", "out_rev", "out_irrev", "out_irrev_fwd", "out_irrev_bwd", "out_closed", "outLP_total", "outLP_rev", "outLP_irrev", "outLP_irrev_fwd", "outLP_irrev_bwd", "outLP_closed"]], B["input_rev"], B["input_bwd"],
                    B["fixIrr_bwd_to_fwd"], B["rev_to_fwd_flipped"], B["rev_to_fwd_bounds"], B["rev_flipped_still_rev"],
                    B["reorient_other"], B["reoriented_at_LP"],
                    B["guard_kept_vs_FVA"], D.get("time_s", ""), A.get("time_s", ""), n,
                    "; ".join(sorted({x["checks"] for x in rows if x["model"] == m}))])
with open(out + "_reoriented.csv", "x", newline="") as fh:
    tc = ["model", "mode", "job", "rxn", "category", "input_class", "lp_class", "flipped", "bounds_changed", "lbIn", "ubIn",
          "lbAtLP", "ubAtLP", "fate", "consistent"]
    w = csv.DictWriter(fh, tc); w.writeheader()
    for m in sorted(G, key=lambda m: (typ(m), m)):
        B = G[m].get("def") or G[m].get("auto")
        for t in B["_TR"]: w.writerow(dict(t, model=m, mode=B["mode"], job=B["job"]))
bad = [r for r in rows if r["checks"] != "OK"]
print(f"{len(rows)} runs, {len(G)} models; runs with a failed check: {len(bad)}")
for r in bad: print("  ", r["model"], r["mode"], r["job"], r["checks"])
print("modes:", dict(collections.Counter(r["mode"] for r in rows)), "| types:", dict(collections.Counter(typ(m) for m in G)))
print("written", out + "_runs.csv", out + "_models.csv", out + "_reoriented.csv")
