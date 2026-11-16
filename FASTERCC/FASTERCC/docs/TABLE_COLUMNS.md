# FasterCC tables: what each column means

16 Nov 2026

## How to read the tables

Three tables describe the 58 real models one row each: Table MP (models and FasterCC preprocessing), Table PM (results per model) and Table 3 extended (FVA reversibility). This document gives the meaning of every column, in the order the columns appear in each file.

**Direction of a reaction, from its bounds** (lb = lower bound, ub = upper bound). Used for every "reversible", "forward", "backward" and "closed" count:

| Class | Rule on the bounds | Meaning |
| --- | --- | --- |
| reversible | lb < 0 < ub | can be given flux in both directions |
| irreversible forward | lb ≥ 0 and ub > 0 | forward only |
| irreversible backward | lb < 0 and ub ≤ 0 | backward only (written "the wrong way round") |
| closed | lb = ub = 0 | cannot carry flux by its bounds |

**FVA terms.** FVA (flux variability analysis) computes, for each reaction, its smallest and largest possible flux in the input model. *FVA-consistent* = FVA says the reaction can carry flux; *FVA-blocked* = it cannot. In Table 3 extended a reaction is *reversible by FVA* when its minimum is below −1e-6 and its maximum above 1e-6.

**Step record.** A FasterCC run that wrote report files saying what happened to every reaction (removed as a dead end, by Poolman, in the LP stage, or kept). The step records used here come from FasterCC build EXA2 (jobs 16199964–16200159). They keep the same reactions as the final benchmark in all 54 models where both exist.

**Special cell values.**

| Value | Meaning |
| --- | --- |
| not recorded | no step record for this model in the files available (Harvey_1_03c, Harvetta_1_03c, Microbiota, ME) |
| no result / run not finished | the method's run did not finish (FasterCC on ME) |
| FVA file not available | the FVA minimum and maximum per reaction are on the HPC but not yet copied here |
| – (in the .md and HTML views) | not available or not applicable |
| empty cell (CSV) | the column does not apply to this row, for example final-benchmark columns for the cluster networks |

## Table MP: models and FasterCC preprocessing

File `table_models_preprocessing.csv` (also .xlsx, .md). One row per real model (58) and per cluster network with no reaction closed (Cluster01–09). The columns follow one reaction set from the input file, through FasterCC's steps, to the final output.

**The model**

| Column | Meaning |
| --- | --- |
| model | model name, as in the benchmark |
| type | model group: Small test, AGORA microbes, Curated microbial, Human generic, Fetal tissues, Whole-body, Microbiota, ME model, or Cluster network |
| reactions | number of reactions in the input model (columns of the stoichiometric matrix S) |
| metabolites | number of metabolites (rows of S), as written in FasterCC's run log; "not recorded" when there is no step record |
| epsilon | threshold ε used in the benchmark: the smallest flux a reaction must reach to count as able to carry flux (1e-4 for most models, 1e-6 for the whole-body models, 1 for the ME model, as set in the task list) |

**The input bounds** (classes defined in the first section)

| Column | Meaning |
| --- | --- |
| input_reversible | reactions with lb < 0 < ub in the input |
| input_irreversible_forward | reactions with lb ≥ 0 and ub > 0 |
| input_irreversible_backward | reactions with lb < 0 and ub ≤ 0 |
| input_closed | reactions with lb = ub = 0 |

**What FasterCC removes, step by step** (the four removal columns plus kept_consistent add up to `reactions`; the script stops if they do not)

| Column | Meaning |
| --- | --- |
| removed_dead_ends | reactions removed as dead ends in the first preprocessing stages (P1–P2): a reaction is a dead end when one of its metabolites is used by only that reaction (type I), or when all reactions of one of its metabolites are irreversible in the same direction, so the metabolite can only be made or only be used (type II) |
| removed_poolman_blocked | reactions the exact Poolman test proves blocked (stage P3): the test looks for a reaction whose row is zero in every null-space basis of S, computed in exact arithmetic by the Java engine |
| removed_dead_ends_after_poolman | reactions that become dead ends only after the Poolman-blocked reactions are removed (stage P3) |
| removed_lp_stage | reactions removed in the LP stage: by the batch LP7, by a singleton LP3, skipped by an exact certificate after the batch, or removed later by the guarded closure (sum of these four fates, automatic mode) |
| kept_consistent | reactions FasterCC keeps as consistent (able to carry flux at ε); for the four models with no step record, the count from the final benchmark |

**How the preprocessing reorients reactions** (for the network handed to the LP stage)

| Column | Meaning |
| --- | --- |
| reversible_made_irreversible | input-reversible reactions that the preprocessing made irreversible because they can carry flux in one direction only; = the next two columns added |
| of_which_backward_only_flipped | … of these, the ones that can only run backward: their column of S is flipped (multiplied by −1) so they become forward-irreversible |
| of_which_forward_only_bounds | … of these, the ones that can only run forward: their lower bound is set to 0, no flip |
| backward_irreversible_written_forward | input reactions that were already backward-irreversible (lb < 0, ub ≤ 0) and are written forward by flipping their column (FasterCC's fixIrr step) |
| reversible_flipped_still_reversible | reversible reactions whose column was flipped but that stay reversible |
| columns_flipped_total | all reactions whose column of S was flipped, of any kind |

**Reversible reactions at three points**

| Column | Meaning |
| --- | --- |
| reversible_at_lp | reversible reactions in the network handed to the LP stage (after dead-end and Poolman removal and after reorientation) |
| reversible_in_output_step_record | reversible reactions among the kept reactions, with FasterCC's reoriented bounds, in the step-record run |
| reversible_in_output_final | the same count in the final benchmark's FasterCC output model (its returned bounds) |

**Checks against the final benchmark**

| Column | Meaning |
| --- | --- |
| kept_final_benchmark | reactions FasterCC kept in the final benchmark run |
| kept_same_as_final | "yes" when the step record keeps the same number of reactions as the final benchmark (true for all 54 models with both) |
| fva_consistent | reactions the benchmark's FVA check of the input model says can carry flux |
| fasterCC_kept_fva_blocked_final | reactions FasterCC kept in the final benchmark although FVA calls them blocked (mostly reactions kept by FasterCC's guard, for example 73 in each fetal tissue) |

**Provenance**

| Column | Meaning |
| --- | --- |
| input_conversion | what was changed in the input file before any method ran, in plain words (for example "missing constraint senses set to equality", or "A matrix split into S and coupling constraints C"); "none" when the model was already in the form S·v = 0 |
| step_record | which run the step columns come from (build, mode, job number), or "none in the files available" |

## Table PM: results per model

File `table_per_model.csv` (also .xlsx, .md). One row per real model (58), from the final benchmark. The columns that end in `_FasterCC`, `_FASTCC` and `_SPECTRA` give the same quantity for each method.

| Column | Meaning |
| --- | --- |
| model | model name |
| type | model group (as in Table MP) |
| reactions | number of reactions in the input model |
| fva_consistent | reactions the benchmark's FVA check of the input model says can carry flux: the reference the methods are compared with |
| kept_FasterCC, kept_FASTCC, kept_SPECTRA | reactions in the consistent network each method returned |
| removed_consistent_… | reactions FVA says can carry flux but the method removed. This is an error: the network loses a working reaction |
| kept_blocked_… | reactions FVA says cannot carry flux but the method kept. For FasterCC these are mostly reactions kept on purpose by its guard |
| time_s_… | run time in seconds from the method's own timer (results_real). The Slurm elapsed times in Table R6 differ by a few seconds to minutes. "no result" = run did not finish |
| reversible_pct_FasterCC | reversible reactions as % of the reactions FasterCC kept, with FasterCC's returned (reoriented) bounds |
| reversible_pct_FASTCC, reversible_pct_SPECTRA | the same %, with the input bounds (these methods do not reorient) |
| notes | anything unusual in the comparison, for example output reaction IDs that are not in the input model |

## Table 3 extended: FVA reversibility

File `table_fva_reversibility.csv` (also .xlsx, .md). One row per real model, per cluster network with no reaction closed, and for iMM1415 and iMM1865 (preprint Table 3). Only Recon3D_mat, HumanGEM12971, iMM1415 and iMM1865 have FVA columns so far; the others fill in when their FVA files are copied from the HPC.

| Column | Meaning |
| --- | --- |
| model, type, reactions | as in Table MP |
| fva_consistent | reactions FVA says can carry flux (from the final benchmark) |
| fva_file | the FVA file used, or "FVA file not available" |
| fva_reversible | reactions FVA says can carry flux in both directions (minimum < −1e-6 and maximum > 1e-6) |
| fva_forward_only | can carry flux forward only (maximum > 1e-6, minimum ≥ −1e-6) |
| fva_backward_only | can carry flux backward only (minimum < −1e-6, maximum ≤ 1e-6) |
| fva_blocked | cannot carry flux in either direction |
| reversible_by_bounds | reversible reactions by the input bounds (lb < 0 < ub) |
| made_irreversible_by_preprocessing | input-reversible reactions FasterCC's preprocessing made irreversible (same as reversible_made_irreversible in Table MP, default-mode run) |
| same_direction_in_fva | … of these, the ones FVA agrees with, checked by reaction ID: a reaction made forward has no backward flux in FVA, one flipped has no forward flux |
| other_direction_in_fva | … the ones FVA disagrees with (FasterCC closed a direction FVA allows). 0 is correct |
| ids_not_in_fva | reoriented reactions whose ID is not in the FVA file (0 means every one was checked) |
| reversible_in_fasterCC_network | reversible reactions left in FasterCC's consistent network after preprocessing (kept reactions, reoriented bounds) |
| gap_left_reversible_but_one_direction_in_fva | reversible_in_fasterCC_network minus fva_reversible: reactions FasterCC leaves reversible although FVA says they can run in one direction only. This is what FasterCC's optional reversibility check (revopt = 1) is meant to close |
| reversible_in_final_benchmark_output | reversible reactions in the final benchmark's FasterCC output (as reversible_in_output_final in Table MP) |
| optional_reversibility_check | whether the optional check was run: "not run with the final version" for every model so far |

The gap is a difference of counts, not a list of reaction IDs. It equals the number of one-direction reactions left reversible only because no reaction FasterCC made irreversible goes against FVA (other_direction_in_fva = 0).

## Where each column comes from

All files are in folder `06_method_comparison_tables` of `FasterCC_results_reproducible_16Nov.zip`. Each table is rebuilt by its script, which refuses to overwrite an existing file.

| Table | Script | Input files |
| --- | --- | --- |
| MP (models and preprocessing) | `make_table_models_preprocessing.py` | `fcc_steps_EXA2_16Nov_models.csv` and `_runs.csv` (step records, read from each run's REPORT files by `fcc_step_table.py`), `methods_real_16Nov_models.csv` (final benchmark), `real_runs_16Nov.tsv` (input conversion), `tasks_real.tsv` (ε) |
| PM (results per model) | `make_table_per_model.py` | `methods_real_16Nov_models.csv` (counted by `compare_methods.py` from the TSV files written by `export_compare.m`) |
| 3 extended (FVA reversibility) | `make_table_fva_reversibility.py` | `fcc_steps_EXA2_16Nov_models.csv`, `fcc_steps_EXA2_16Nov_reoriented.csv` (every reaction reoriented, default mode), `methods_real_16Nov_models.csv`, `fva_lbub/` (FVA minimum and maximum per reaction) |

The model type in all three tables comes from `02_FasterCC_vs_SPECTRA/out/figSP3_real_model_times.csv`.
