# FASTCORE, correctness, time and Table 3: figure pack, 16 Nov 2026

Rebuild: `python3 -I make_core_correct_figs.py data out` (Python 3, matplotlib, openpyxl; never overwrites).

| Item | Rebuilt as | Source data |
|---|---|---|
| Preprint Figure 4 / ledger Figure R3 | figP4_R3_fastcore_preprocessing | CORE_by_task.csv (250 FASTCORE tasks, clusters 1-10, CPLEX) |
| Preprint Table 2 | table2_fastcore_top5 | CORE_by_task.csv |
| Preprint Supplementary Table 2 | tableS2_fastcore_by_size_and_core | CORE_by_task.csv |
| Ledger Table R7 | tableR7_fastcore_by_cluster | CORE_by_task.csv |
| Ledger Figure R5 | figR5_correctness | methods_real_16Nov_models.csv, methods_bench_tasks.csv, fva_check/*.json |
| Ledger Table R6 | tableR6_wall_and_cpu_time | sacct_largest.txt (Slurm), tasks_real.tsv |
| Ledger Table R22 | tableR22_objective | BIGOBJ_16226422.out, BIGOBJ2_16210166.out |
| Preprint Table 3 | table3_reversible_reactions, table3_check_by_reaction_id | data/table3 (FVA files; FasterCC 6 Nov runs) |

Notes:
- FASTCORE benchmark = 250 tasks of the final run (jobs 16207705-09, FasterCC AW2, CPLEX); 10 clusters x 5 core sizes x 5
  replicates by design. The expanded networks start at 12,160 reactions (preprint: 10,600).
- Table 3: bounds and FVA columns reproduce the preprint exactly. The FasterCC column here is preprocessing only (no optional
  reversibility check); the check will be run after the 16 Nov maintenance and added.
- Figure R5b compares counts: the FVA consistency check stores counts, not IDs; it covers clusters 1-6 so far.
Colours (Okabe-Ito): FasterCC #0072B2, FASTCC #E69F00, SPECTRA #009E73.
