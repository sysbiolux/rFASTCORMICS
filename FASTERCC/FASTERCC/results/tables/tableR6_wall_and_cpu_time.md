**Table R6 | Elapsed and CPU time on the largest real models (Slurm accounting).** Elapsed and TotalCPU of each run's Slurm
job (sacct), cores used = CPU / elapsed. Each method ran in its own job with the same allocation per model class (38 CPUs
for the whole-body, Microbiota and ME models, 10 for the human models). The ME FasterCC run completed in Slurm but has no
exported result (ME restart), so its time is listed but not compared elsewhere. Sources: sacct_largest.txt, tasks_real.tsv.

| model | CPUs allocated | elapsed s FasterCC | elapsed s FASTCC | elapsed s SPECTRA | CPU s FasterCC | CPU s FASTCC | CPU s SPECTRA | cores used FasterCC / FASTCC / SPECTRA | FASTCC / FasterCC: elapsed -> CPU | SPECTRA / FasterCC: elapsed -> CPU | Slurm jobs |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Microbiota | 38 | 9,123 | 342,680 | 1,446 | 28,856 | 1,838,619 | 4,522 | 3.2 / 5.4 / 3.1 | 37.6 -> 63.7 | 0.16 -> 0.16 | 16207712_15 16172902_13 16172902_14 |
| Harvetta_aafix_E6 | 38 | 3,211 | 118,801 | 396 | 19,217 | 2,071,486 | 4,078 | 6.0 / 17.4 / 10.3 | 37.0 -> 107.8 | 0.12 -> 0.21 | 16207712_12 16172902_10 16172902_11 |
| Harvey_aafix_E6 | 38 | 2,821 | 115,865 | 393 | 35,435 | 2,039,614 | 4,554 | 12.6 / 17.6 / 11.6 | 41.1 -> 57.6 | 0.14 -> 0.13 | 16207712_9 16172902_7 16172902_8 |
| Harvey_1_03c | 38 | 83,262 | 111,275 | 325 | 82,985 | 2,005,454 | 3,729 | 1.0 / 18.0 / 11.5 | 1.3 -> 24.2 | 0.00 -> 0.04 | 16207712_3 16172902_1 16172902_2 |
| Harvetta_1_03c | 38 | 85,951 | 93,452 | 475 | 85,614 | 1,604,895 | 5,076 | 1.0 / 17.2 / 10.7 | 1.1 -> 18.7 | 0.01 -> 0.06 | 16207712_6 16172902_4 16172902_5 |
| ME_dbl_E0_L1_00e_2 | 38 | 193,736 | 14,746 | 249 | 265,093 | 227,116 | 3,114 | 1.4 / 15.4 / 12.5 | 0.1 -> 0.9 | 0.00 -> 0.01 | 16209069_18 16172902_16 16172902_17 |
| Recon3D | 10 | 34 | 1,540 | 41 | 117 | 10,260 | 75 | 3.4 / 6.7 / 1.8 | 45.3 -> 87.7 | 1.21 -> 0.64 | 16207711_30 16172901_28 16172901_29 |
| Human1 | 10 | 49 | 1,514 | 51 | 168 | 10,566 | 222 | 3.4 / 7.0 / 4.4 | 30.9 -> 63.0 | 1.04 -> 1.32 | 16207711_39 16172901_37 16172901_38 |
| HumanGEM | 10 | 39 | 1,180 | 47 | 160 | 8,554 | 197 | 4.1 / 7.2 / 4.2 | 30.3 -> 53.5 | 1.21 -> 1.23 | 16207711_33 16172901_31 16172901_32 |
