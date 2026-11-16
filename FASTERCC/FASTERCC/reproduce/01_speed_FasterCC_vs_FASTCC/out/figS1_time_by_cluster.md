**Figure S1 | Run time of FasterCC, FASTCC and SPECTRA on the cluster benchmark.** Median wall time per consistency-check
task (60 tasks per cluster: removal levels 0, 10, 50, 100, 200 and 500 reactions, 10 replicates each) for the 19 cluster
networks (21,200 to 183,920 reactions), with CPLEX (left) and Gurobi (right); log scale. Every method ran on the same input
files; with Gurobi all three ran in the same job, with CPLEX FasterCC ran in one job and FASTCC and SPECTRA in another;
the 15 SPECTRA tasks lost to a licence failure are the SPEC15 reruns. Over the 1,140 tasks the geometric mean of the per-task speed-up FASTCC/FasterCC is 8.11x
with CPLEX and 4.51x with Gurobi; SPECTRA/FasterCC is 0.54x and 1.33x.
Data: figS1_time_by_cluster.csv; raw times: data/speed_tasks.tsv; code: make_speed_figs.py.
