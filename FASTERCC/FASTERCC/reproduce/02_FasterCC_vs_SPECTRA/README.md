# FasterCC compared with SPECTRA (and FASTCC): figure pack, 16 Nov 2026

Rebuild: `python3 -I make_spectra_figs.py data out` (Python 3, matplotlib; writes nothing over existing files).
Each figure: `.svg`, `.pdf`, `.png` (300 dpi), plotted data `.csv`, caption `.md`. Each table: `.csv` + `.md`.

| Item | What it shows | Data |
|---|---|---|
| figSP1 | time per task, FasterCC vs SPECTRA, 1,140 cluster tasks x 2 solvers | speed_tasks.tsv |
| figSP2 | median SPECTRA/FasterCC time ratio, every network x number of blocked reactions | speed_tasks.tsv |
| figSP3, tableSP2 | run times of the three methods on the 58 real models | methods_real_16Nov_models.csv |
| figSP4 | real models: reactions removed although FVA-consistent / kept although FVA-blocked | methods_real_16Nov_models.csv |
| figSP5 | reversible reactions left in the output (FasterCC after reorientation vs input directions) | methods_real_16Nov_models.csv, methods_bench_tasks.csv |
| tableSP1 | number of reactions kept, method vs method, per task, after the SPEC15 reruns | methods_bench_tasks.csv, spec15_tasks.tsv |

Sources: `speed_tasks.tsv` (attached 12 Nov; SPECTRA times of the 15 licence-failed tasks are the SPEC15 reruns);
`methods_*` from `export_compare.m` v2 (md5 47aaf1bd…) + `compare_methods.py` (md5 83d80a61…) on the cluster result
files; `spec15_tasks.tsv` copied from the cluster output of 16 Nov 03:46 (`compare_spec15/bench_tasks.tsv`).
Real-model times are single runs as stored in `results_real`; the solver used for them is not a column of the export.
Colours (Okabe-Ito): FasterCC #0072B2, FASTCC #E69F00, SPECTRA #009E73.
