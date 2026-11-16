**Figure SP3 | Run time of FasterCC, FASTCC and SPECTRA on the 58 real models.** Models grouped (italic labels) and ordered by
size within each group; reactions in brackets. FasterCC has no time for ME_dbl_E0_L1_00e_2 (its run is the ME restart still
running). Times as stored in results_real (one run per model and method). On small and mid-size real models SPECTRA is
usually the fastest; FasterCC is faster than FASTCC on the larger models. The FasterCC times of Harvey_1_03c and Harvetta_1_03c (about 23 h) are dominated by careful-mode re-solves: these two models are infeasible with their own bounds at 1e-9, so about 850 LPs fall back to single-thread primal simplex (tableR6b_careful_mode_1_03c); their feasible aafix versions take 47-54 min. Per-group geometric means: tableSP2_real_times_by_group.
Data: figSP3_real_model_times.csv; code: make_spectra_figs.py.
