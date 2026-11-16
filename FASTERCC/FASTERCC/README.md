# FasterCC: benchmark results (November 2026)

FasterCC returns the flux-consistent part of a constraint-based metabolic model: the reactions that can carry a flux of at
least ε at steady state. Before any linear program (LP), it removes reactions that cannot carry flux (dead ends, and
reactions blocked by the stoichiometry alone, proven by an exact null-space analysis test (Poolman), then solves far fewer LPs than FASTCC.

This folder holds the results of the November 2026 benchmark of the new FasterCC against FASTCC and SPECTRA:
- 1,140 cluster tasks (CPLEX and Gurobi), 250 FASTCORE tasks and 58 real models;
- every table and figure, with the data and scripts that rebuild them;
- the consistent network each method returned for each real model (as reaction lists);
- two toy models.

**The code.** The FasterCC version of the preprint is in the parent folder (`FASTERCC/`). The code of the new version will
be added here once its cleaned-up release has passed the regression tests against the benchmarked version.

---

## Contents

| Folder | What is in it |
|---|---|
| `results/tables/` | every results table (CSV, Markdown; the three per-model tables also as Excel) |
| `results/figures/` | every figure (PNG) |
| `results/consistent_networks/` | for each of the 58 real models, which reactions FasterCC, FASTCC and SPECTRA kept, and FasterCC's output bounds (one TSV per model; see its README) |
| `reproduce/` | data and scripts that rebuild every table and figure, one folder per figure type |
| `docs/` | `FIGURE_GUIDE.md` (description and message of every figure) and `TABLE_COLUMNS.md` (meaning of every table column) |
| `toy_examples/` | the two toy models, MATLAB scripts that run the methods on them, and a Python check needing no solver licence |

---

## Benchmark results (8–11 Nov 2026)

Every number below is counted from the files in `reproduce/`. Each table and figure has a description and a message in `docs/FIGURE_GUIDE.md`; every table column is explained in `docs/TABLE_COLUMNS.md`.

**Same answer as FASTCC.** On the 1,140 cluster tasks (19 networks of 21,200–183,920 reactions, 0–500 reactions closed at random, 10 replicates), FasterCC keeps the same number of reactions as FASTCC:
- in 1,140 of 1,140 tasks with Gurobi;
- in 1,138 with CPLEX.

Every reaction-ID difference between the methods is one biomass reaction, decided by exact LPs (Table R23).

**Faster than FASTCC.** Per-task speed-up, FASTCC time ÷ FasterCC time on the same network (geometric mean):

| | CPLEX | Gurobi |
|---|---|---|
| all 1,140 tasks | 8.11× | 4.51× |
| no reaction closed | 4.0× | 2.7× |
| 500 reactions closed | 25.5× (per task 10.4–59.7×) | 11.2× (4.6–56.4×) |

On large real models FasterCC is 37–45× faster than FASTCC in elapsed time:

| Model | FasterCC | FASTCC |
|---|---|---|
| Harvey_aafix (whole body) | 2,821 s | 115,865 s |
| Harvetta_aafix (whole body) | 3,211 s | 118,801 s |
| Microbiota (173,353 reactions) | 9,123 s | 342,680 s |
| Recon3D | 34 s | 1,540 s |

**Against SPECTRA.** With CPLEX SPECTRA is about twice as fast (median ratio 0.54–0.58). With Gurobi FasterCC is faster on every network from 66,400 reactions, and in 822 of 1,140 tasks.

SPECTRA's speed comes with errors:
- it keeps FVA-blocked reactions on 7 large models, e.g. Microbiota 954 and ME 2,453;
- on Microbiota only FasterCC's network can still produce community biomass;
- its output can change when SPECTRA is run on it again.

**Correct on the real models.** Against FVA:
- **FasterCC** removes a reaction that can carry flux on only 2 models: the two Harvetta models, 3 reactions each, the same 3 that FASTCC removes.
- **FASTCC** does so on 10 models.
- **Blocked biomass reactions kept:** FasterCC 1 (CPLEX) and 0 (Gurobi), FASTCC 2 and 0, SPECTRA 15 and 2.
- **No method** removes a working biomass reaction.

**Biomass check** (`reproduce/04_biomass_wholebody_ME/`, details in its README). Exact LPs decided every biomass reaction of every cluster task's input network: 37,620 rows, 37,523 that can carry flux and 97 blocked. Each verdict was then compared with what each method did with that reaction.
- 58 of the 97 blocked ones are closed by their own bounds, and every method removes them.
- 39 have open bounds but no route left to feed them. Every wrongly kept biomass reaction, for every method, is one of these 39.
- FasterCC keeps one of them, in task CC_c02_l050_r09, with CPLEX only. Every CPLEX method keeps that one.
- FASTCC keeps that one and one more. SPECTRA keeps 15 with CPLEX and 2 with Gurobi.

`python3 -I check_biomass.py` in that folder recounts Table R23 from the raw rows and lists every case.

**Better input for downstream methods.**
- FasterCC's output has about half as many reversible reactions as FASTCC's and SPECTRA's.
- FASTCORE is faster on the network cleaned by FasterCC in 211 of 250 tasks.

**Known limits.**
- On inputs that are infeasible with their own bounds at 1e-9 (the original whole-body models `_1_03c`), FasterCC's precision safeguard ("careful mode") takes about 23 h (Table R6b).
- FasterCC's run on the ME model has not finished.

### Tables (`results/tables/`)

| File | Content |
|---|---|
| `table_per_model.*` | the 58 real models: kept, errors against FVA, time and reversible share, per method |
| `table_models_preprocessing.*` | each model and what FasterCC removes at each step (dead ends, exact Poolman, LP stage) |
| `table_fva_reversibility.*` | Table 3 extended: reversible reactions by FVA against FasterCC, and the gap left for the optional check |
| `table1_*`, `tableS1_*`, `table_speed_by_cluster.*` | speed against FASTCC (preprint Table 1 and Supplementary Table 1 layout) |
| `tableSP1_*`, `tableSP2_*` | agreement and run times against SPECTRA |
| `tableR23_biomass.*` | blocked and working biomass reactions kept or removed by each method (exact LPs) |
| `tableR6_*`, `tableR6b_*` | elapsed and CPU time on the large real models; careful mode on `_1_03c` |
| `tableR22_objective.*` | whether the reduced networks still reach their objective |
| `tableWB_spectra_wholebody.*` | SPECTRA on the whole-body models at ε 1e-4 and 1e-6 |
| `table2_*`, `tableR7_*`, `tableS2_*` | FASTCORE with and without FasterCC |
| `table3_*` | preprint Table 3: reversible reactions by bounds, after FasterCC and by FVA |
| `tableME_snapshot.*` | the ME model: where FasterCC's unfinished run stands |

---

## Rebuilding the tables and figures (`reproduce/`)

Each folder holds its data (`data/`), its script and its outputs (`out/`). The scripts need Python 3 with numpy, matplotlib and openpyxl, and they never overwrite an existing output folder.

```bash
cd reproduce/01_speed_FasterCC_vs_FASTCC         && python3 -I make_preprint_figs.py data new_out
cd ../02_FasterCC_vs_SPECTRA                     && python3 -I make_spectra_figs.py data new_out
cd ../03_FASTCORE_correctness_time_objective_Table3 && python3 -I make_core_correct_figs.py data new_out
cd ../04_biomass_wholebody_ME                    && python3 -I make_bio_wb_me.py data new_out
cd ../05_loopless                                && python3 -I plot_loopless.py data/ll_15Nov_permodel.csv new_out
cd ../06_method_comparison_tables                && python3 -I make_table_per_model.py new_per_model.csv
```

`MD5SUMS.txt` lists the checksum of every file in this folder (`md5sum -c MD5SUMS.txt`, run inside the folder).

---

## Toy examples (`toy_examples/`)

Two small networks show what each method gets right or wrong, and why.

You can run each toy without MATLAB or a commercial solver:

```bash
python3 toy_examples/check_toys_python.py      # needs numpy and scipy
```

The script prints FVA, the exact Poolman result, and the consistent set at ε, closed by repeating until nothing changes. Its output is saved in `check_toys_python_output.txt`.

### Toy 1: a reaction-level ε threshold keeps a producer but removes its consumer (`toy_spectra_fastcc.m`)

The toy model has two branches from A[c]:

```
v1:  -> A[e]                 (upper bound 1e-4 = epsilon)
v2:  A[e] -> A[c]
B branch:  v3: A[c] -> B[c]      v4: B[c] -> 0.1 B[e]      v5: B[e] ->
C branch:  v6: A[c] -> C[c]      v7: C[c] -> C[e]          v8: C[e] ->
```

**Why v5 is never selected.** Because v1.ub = 1e-4 and ε = 1e-4, the maximum flux entering the model is 1e-4. In the B branch, even if v4 carries 1e-4, the maximum possible flux through v5 is only 0.1 × 1e-4 = 1e-5. That is below ε, so v5 is not selected as consistent.

**What SPECTRA does** (`spectraCC(model, 1e-4, 'stoichiometry')`, from our run of the script in this folder):

| Run | Reactions kept |
|---|---|
| First run | v1, v2, v3, v4, v6, v7, v8 |
| Run again on its own output | v1, v2, v6, v7, v8 |

The first run removed v5 but kept v3 and v4. In that output v4 produces B[e] and nothing consumes it, so B[e] is a dead end and the B branch is inconsistent. Run again, SPECTRA removes v3 and v4. A consistent output must not change when the method is run on it again; SPECTRA's first answer does.

**The correct answer at ε.** It is {v1, v2, v6, v7, v8}. The Python check reaches it independently:
1. FVA: v5 can carry at most 1e-5.
2. Remove v5; FVA then gives v3 and v4 no flux.
3. Remove them; nothing else changes.

`toy_spectra_fastcc.m` prints, for SPECTRA, FASTCC and FasterCC:
- the set each keeps;
- the dead-end metabolites left in its output;
- whether running the method again on its own output changes the set.

**The same mechanism in a real model.** We traced this behaviour in the Harvey whole-body model, in the dolichol pathway. Reactions with different stoichiometric scaling carry fluxes of different magnitudes. Upstream producer reactions pass the ε threshold, while the downstream consumers that balance them stay below ε and are removed. This leaves dead ends in SPECTRA's "consistent" model. Dead ends and structural consistency were analysed with FasterCC's own functions `findType1and2deadEnd.m`, `structureAnalyseFASTERCC.m` and `fixIrrFASTERCC.m`. On the whole-body models at ε = 1e-6, SPECTRA keeps 62–70 reactions that FVA calls blocked, or 430–558 with the coupling constraints (`results/tables/tableWB_spectra_wholebody.*`). It also removes further reactions when run again on its own output.

### Toy 2: what exact Poolman removes, and what it must keep (`toy_poolman.m`)

```
R1: -> A      R2: A -> B      R3: B ->        pathway: carries flux
R4: A -> C    R5: C -> D      R6: D -> C      R4 feeds the pool C + D, which nothing drains
R7: B -> E                                    E is never consumed (dead end)
all irreversible, bounds 0..10
```

- **R4 is blocked by the stoichiometry alone.** Adding the balances of C and D gives d(C + D)/dt = v4, so v4 = 0 at steady state. Its row is zero in every null-space basis of S, and y = C + D is its exact certificate (Sᵀy = e_R4).
- **R7 is blocked:** E is a dead end.
- **R5 and R6 can carry flux.** They form a free loop (v5 = v6 balances both C and D), so FASTCC and FVA call them consistent, and FasterCC must keep them.

FasterCC's exact Poolman removes R4 and R7 before any LP is solved; the Python check gives the same two. The script also runs FVA and the methods it finds on the MATLAB path.
