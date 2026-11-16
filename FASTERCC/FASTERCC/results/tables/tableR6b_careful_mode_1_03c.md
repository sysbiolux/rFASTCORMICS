**Why FasterCC took 23 h on the `_1_03c` whole-body models (from the run logs, BENCH_real_NEW_large 16207712_3, _6, _9; the `[COUNT]` timers are inclusive, so lines overlap and do not add up).** The exact Poolman step is the same on all three runs (about 2 s per call, 3 calls, about 1,000 reactions removed). All three runs switch FasterCC's precision safeguard ("careful mode") on (`FCC_PRECISION_MODE=auto`: 22.3 digits needed, threshold 18): each LP is first tried at a strict tolerance (1e-9) and the answer checked against all constraints. The `_1_03c` networks are infeasible with their own bounds at 1e-9 (BIGOBJ), so the strict set-ups find no solution and the LP falls back to set-up 3 (primal simplex at 1e-6), which CPLEX runs on one thread.

| | Harvey_1_03c | Harvetta_1_03c | Harvey_aafix_E6 |
|---|---|---|---|
| FasterCC total | 83,233 s | 85,914 s | 2,788 s |
| Strict set-ups without solution (dual simplex 1e-9 / barrier 1e-9) | 848 / 848 | 850 / 850 | 0 / 0 |
| Answers above CPLEX's feasibility tolerance when checked | 639 | 638 | 0 |
| LPs answered by set-up 3 (primal 1e-6): calls, time | 848 LPs, 81,924 s (97 s each) | about 850 (similar) | – (172 LPs by set-up 1, barrier 1e-9: 2,091 s, 12 s each) |
| LP3 (all calls) | 495 calls, 65,582 s (132 s each) | – | 1,182 calls, 845 s (0.7 s each) |
| LP7 batch (all calls) | 144 calls, 13,605 s (94 s each) | – | 92 calls, 1,404 s (15 s each) |
| Careful-mode proofs | 105 calls, 4,245 s | – | 45 calls, 179 s |
| Cores used (Slurm CPU / elapsed) | 1.0 | 1.0 | 12.6 |
| Result | 68,233 kept; FVA: 0 consistent removed | 70,355 kept; 3 removed (same 3 as FASTCC) | 68,072 kept; 0 removed |

For the paper: the fair speed comparison on whole-body models is the `aafix` versions (feasible): FasterCC 37–41× faster than FASTCC in elapsed time (58–108× in CPU time). The `_1_03c` versions are reported as the case where the model is numerically infeasible with its own bounds: careful mode re-solves about 850 LPs on one thread, so FasterCC is only 1.1–1.3× faster than FASTCC in elapsed time (19–24× in CPU time), with a correct result. Planned after the maintenance: one run of Harvey_1_03c with careful mode off, to measure its cost directly and compare the kept set.
