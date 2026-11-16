**Figure 4 (preprint) / R3 (ledger) | FASTCORE with and without FasterCC preprocessing.** a, Median (line) and interquartile
range (band) of FASTCORE runtime on the raw expanded cluster network and on the same network after FasterCC cleaning, by
core size (50 tasks per core size: clusters 1-10 of 12,160 to 93,520 reactions x 5 replicates; CPLEX); dotted grey, cleaning
plus one FASTCORE run. b, Per-task speed-up (time without / with preprocessing), median and IQR; dotted line, no gain.
FASTCORE was faster on the cleaned network in 211 of 250 tasks and at least 10% faster (time ratio >= 1/0.9) in 182.
The core was kept in all 250 tasks with both inputs; the two outputs were never identical (Jaccard index in
tableR7_fastcore_by_cluster). Note: the preprint's networks started at 10,600 reactions; these expanded inputs start at
12,160 (cluster 1). Data: figP4_R3_fastcore_preprocessing.csv; source CORE_by_task.csv; code: make_core_correct_figs.py.
