#!/usr/bin/env python3
"""check_toys_python.py -- independent reference for the two toy models, without MATLAB, COBRA or a commercial solver.
For each toy it prints:
  1. FVA: the smallest and largest flux of every reaction at steady state (S v = 0, lb <= v <= ub), by scipy's HiGHS LP;
  2. which reactions can carry flux at all (max |v| > 0) and which reach the threshold epsilon (max |v| >= epsilon);
  3. exact Poolman: reactions whose row is zero in every null-space basis of S, computed in exact fractions
     (a reaction is stoichiometrically blocked <=> there is an exact y with S'y = e_i);
  4. the consistent set at epsilon, closed: remove every reaction whose largest flux is below epsilon, recompute FVA on
     what is left, and repeat until nothing changes (the set a correct method must return; running it again changes nothing).
Usage: python3 check_toys_python.py            (needs numpy and scipy)"""
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog

def fva(S, lb, ub):
    n = S.shape[1]; lo, hi = [], []
    for j in range(n):
        c = np.zeros(n)
        for sgn, out in ((1, lo), (-1, hi)):
            c[:] = 0; c[j] = sgn
            r = linprog(c, A_eq=S, b_eq=np.zeros(S.shape[0]), bounds=list(zip(lb, ub)), method="highs")
            out.append(sgn * r.fun if r.status == 0 else float("nan"))
    return np.array(lo), np.array(hi)

def exact_poolman(Sfr):
    """reactions i with an exact y, S'y = e_i  <=>  e_i in the row space of S  <=>  rank([S; e_i]) = rank(S)"""
    def rank(M):
        M = [row[:] for row in M]; r = 0; ncol = len(M[0])
        for c in range(ncol):
            p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
            if p is None: continue
            M[r], M[p] = M[p], M[r]
            for i in range(len(M)):
                if i != r and M[i][c] != 0:
                    f = M[i][c] / M[r][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
            r += 1
        return r
    r0 = rank(Sfr); n = len(Sfr[0]); out = []
    for i in range(n):
        e = [Fr(0)] * n; e[i] = Fr(1)
        if rank(Sfr + [e]) == r0: out.append(i)
    return out

def report(name, mets, rxns, S, lb, ub, eps):
    print("=" * 100); print(name); print("=" * 100)
    lo, hi = fva(np.array(S, dtype=float), lb, ub)
    blocked_exact = exact_poolman([[Fr(x).limit_denominator(10**12) for x in row] for row in S])
    print(f"{'reaction':10s} {'min flux':>12s} {'max flux':>12s} {'can carry flux':>15s} {'reaches eps':>12s} {'exact Poolman':>14s}")
    for j, r in enumerate(rxns):
        m = max(abs(lo[j]), abs(hi[j]))
        print(f"{r:10s} {lo[j]:12.3e} {hi[j]:12.3e} {('yes' if m > 1e-12 else 'NO'):>15s} {('yes' if m >= eps else 'NO'):>12s} "
              f"{('BLOCKED' if j in blocked_exact else '-'):>14s}")
    print(f"epsilon = {eps:g}; 'reaches eps' uses max |v| >= epsilon")
    keep = list(range(len(rxns))); it = 0
    while True:
        it += 1
        Sk = np.array(S, dtype=float)[:, keep]
        lo2, hi2 = fva(Sk, [lb[j] for j in keep], [ub[j] for j in keep])
        drop = [keep[k] for k in range(len(keep)) if max(abs(lo2[k]), abs(hi2[k])) < eps]
        print(f"  closure pass {it}: below epsilon {[rxns[j] for j in drop] or 'none'}")
        if not drop: break
        keep = [j for j in keep if j not in drop]
    print(f"  consistent set at epsilon (closed): {[rxns[j] for j in keep]}\n")

# ---- toy 1: the SPECTRA / FASTCC toy (COBRA script in toy_spectra_fastcc.m). Default bounds 0..1000, v1 <= 1e-4
mets = ["A[e]", "A[c]", "B[c]", "B[e]", "C[c]", "C[e]"]
rxns = ["v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8"]
S = [[1, -1, 0, 0, 0, 0, 0, 0],        # A[e]:  v1 -> A[e] -> v2
     [0, 1, -1, 0, 0, -1, 0, 0],       # A[c]:  v2 -> A[c] -> v3, v6
     [0, 0, 1, -1, 0, 0, 0, 0],        # B[c]:  v3 -> B[c] -> v4
     [0, 0, 0, 0.1, -1, 0, 0, 0],      # B[e]:  v4 -> 0.1 B[e] -> v5
     [0, 0, 0, 0, 0, 1, -1, 0],        # C[c]:  v6 -> C[c] -> v7
     [0, 0, 0, 0, 0, 0, 1, -1]]        # C[e]:  v7 -> C[e] -> v8
lb = [0] * 8; ub = [1e-4] + [1000] * 7
report("Toy 1 (SPECTRA / FASTCC toy): uptake v1 <= 1e-4, B leaves as 0.1 B[e]", mets, rxns, S, lb, ub, 1e-4)

# ---- toy 2: the Poolman toy (toy_poolman.m). 5 metabolites, 7 irreversible reactions, bounds 0..10
mets = ["A", "B", "C", "D", "E"]
rxns = ["R1", "R2", "R3", "R4", "R5", "R6", "R7"]
S = [[1, -1, 0, -1, 0, 0, 0],          # A: R1 in; R2, R4 out
     [0, 1, -1, 0, 0, 0, -1],          # B: R2 in; R3, R7 out
     [0, 0, 0, 1, -1, 1, 0],           # C: R4, R6 in; R5 out
     [0, 0, 0, 0, 1, -1, 0],           # D: R5 in; R6 out
     [0, 0, 0, 0, 0, 0, 1]]            # E: R7 in, nothing out (dead end)
report("Toy 2 (Poolman toy): R4 feeds the closed pool C + D, R7 ends in E, R5-R6 is a free loop",
       mets, rxns, S, [0] * 7, [10] * 7, 1e-4)
