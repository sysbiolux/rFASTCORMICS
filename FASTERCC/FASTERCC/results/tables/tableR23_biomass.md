**Table R23 | Biomass reactions kept or removed wrongly, per method.** Every biomass reaction of every cluster task (37,620 rows:
37,523 that can carry flux, 97 blocked), decided by exact LPs in the task's input network (biomass check of 10 Nov,
fill-in 2). 'FasterCC (preprint)' = the earlier release run on 527 tasks. No method removes a biomass reaction that can
carry flux. The one blocked biomass FasterCC keeps with CPLEX is task CC_c02_l050_r09 (biomassreaction_2, kept by every CPLEX
method; the KEPTCHK case). Source: data/biocheck.

| method | solver | tasks with a result (of 1,140) | blocked biomass kept (rows) | working biomass removed (rows) | blocked biomass removed | working biomass kept |
|---|---|---|---|---|---|---|
| FasterCC | CPLEX | 1,140 | 1 | 0 | 96 | 37,523 |
| FASTCC | CPLEX | 1,140 | 2 | 0 | 95 | 37,523 |
| SPECTRA | CPLEX | 1,140 | 15 | 0 | 82 | 37,523 |
| FasterCC (preprint) | CPLEX | 527 | 1 | 0 | 46 | 15,874 |
| FasterCC | Gurobi | 1,140 | 0 | 0 | 97 | 37,523 |
| FASTCC | Gurobi | 1,140 | 0 | 0 | 97 | 37,523 |
| SPECTRA | Gurobi | 1,140 | 2 | 0 | 95 | 37,523 |
