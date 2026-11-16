# Speed figures and table of the cluster benchmark (16 Nov 2026)

Everything here is rebuilt from the raw per-task times with one command:

    python3 make_speed_figs.py data out        # needs Python 3 and matplotlib

* `data/speed_tasks.tsv` — 1,140 tasks (clusters 2–20 × removal levels 0, 10, 50, 100, 200, 500 × 10 replicates), wall time in seconds of FasterCC (`FC`), FASTCC and SPECTRA with CPLEX (`_cpx`) and Gurobi (`_grb`); the 15 SPECTRA tasks lost to a licence failure are their SPEC15 reruns.
* `data/cluster_sizes.csv` — reactions (reversible / irreversible) of each base cluster network.
* `out/` — for each figure: `.svg`, `.pdf`, `.png` (300 dpi), the plotted numbers (`.csv`) and the caption (`.md`); `table_speed_by_cluster.csv/.md`.

Checksums: `MD5SUMS.txt`.

Preprint and ledger layouts (same data, plus `data/methods_bench_tasks.csv` for the same-result column):

    python3 make_preprint_figs.py data out

writes Figure 3 (preprint), Table 1 and Supplementary Table 1 (CPLEX and Gurobi), and Figures R1, R2/S5/S6, R4/Figure 5 and S3 of the ledger.
