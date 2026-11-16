**Figure R5 | Correctness of FasterCC, SPECTRA and FASTCC.** a, Real models: SPECTRA time against FasterCC time (log-log;
dashed, equal time; the two points near 8e4 s on the x axis are Harvey_1_03c and Harvetta_1_03c, slowed by careful-mode re-solves, see tableR6b_careful_mode_1_03c); marker area grows with the number of FVA-blocked reactions SPECTRA keeps (labelled when >= 100). b,
Cluster tasks with an independent FVA reference (257 tasks: the networks of clusters 2-6 that are in both the speed
benchmark and the FVA consistency check, which so far covers clusters 1-6 and two networks of cluster 7; Gurobi 12, threshold
0.99 eps, every LP optimal): number of tasks where a method keeps more reactions than FVA finds consistent; solid CPLEX,
hatched Gurobi. Tasks where a method keeps fewer than FVA (removing a reaction that can carry flux): 0 in total. This
comparison uses counts (the FVA check stores counts, not reaction IDs); the reaction IDs of the differences between methods are
in methods_bench_16Nov_diff_ids_finished.csv (all biomass reactions). Data: figR5_correctness.csv; code: make_core_correct_figs.py.
