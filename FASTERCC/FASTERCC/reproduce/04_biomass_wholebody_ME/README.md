# 04 – Biomass check, SPECTRA on the whole-body models, the ME model

This folder holds three results:
- the **biomass check** (Table R23, Figure R23);
- SPECTRA on the whole-body models (Table WB, Figure WB);
- the ME model snapshot (Table ME).

This README is mostly about the biomass check.

## Why check the biomass reactions

A biomass reaction is the one reaction a modeller cannot afford to lose or keep wrongly:
- If a method **removes a biomass reaction that can carry flux**, the reduced network can no longer grow. Every later analysis on it is wrong.
- If it **keeps a biomass reaction that cannot carry flux**, the network looks able to grow when it is not.

Every reaction-ID difference between FasterCC, FASTCC and SPECTRA in the 1,140 cluster tasks turned out to be one biomass reaction. So the biomass reactions are where the methods disagree, and the check below decides which method is right.

## What was tested

For every cluster task (19 networks × 6 levels of closed reactions × 10 replicates = 1,140 tasks), each biomass reaction of the task's **input network** was tested:
- the network after the task closed its random reactions;
- before any method ran.

That gives 37,620 rows, one per biomass reaction per task. The rows split evenly into three types: `biomassreaction_N`, `biomassmaintenance_N` and `biomassmaintenancenoTrTr_N`, 12,540 each.

Each reaction got a verdict from exact LPs (biomass check of 10 Nov, fill-in 2):
- **CONSISTENT:** it can carry flux. `maxflux` is its largest flux, between 11.3 and 755 for the 37,523 working ones.
- **BLOCKED:** it cannot carry any flux (`maxflux` = 0). There are 97 such rows.

The verdict was then compared with what each method did with that reaction in the same task.

The exact-LP scripts that wrote the verdicts are **not in this repository**; only their output is. The `lps` column (values 3–6) is written by those scripts, so its exact meaning is defined there.

## Files

| File | Content |
|---|---|
| `data/biocheck/bio_cNN_fill2.csv` | one file per cluster (02–20), one row per biomass reaction per task |
| `data/biocheck/kept_blocked_list.txt` | the blocked biomass reactions kept by a CPLEX run (method, task, reaction, maxflux, lps, status) |
| `make_bio_wb_me.py` | builds Table R23 and Figure R23 (plus the whole-body and ME tables) |
| `check_biomass.py` | recounts Table R23 independently from the CSVs, explains why each blocked reaction is blocked, and lists every wrongly kept one |
| `out/check_biomass_output.txt` | the output of `check_biomass.py` |

Columns of `bio_cNN_fill2.csv`:

| Column | Meaning |
|---|---|
| `task` | task ID: `CC_c<cluster>_l<reactions closed>_r<replicate>` |
| `cluster`, `level`, `rep` | the same three numbers, split |
| `reaction` | biomass reaction ID |
| `lb`, `ub` | its bounds in the task's input network (0, 0 = closed by the task) |
| `verdict` | CONSISTENT or BLOCKED, from the exact LPs |
| `maxflux` | its largest flux in the exact LP (0 when blocked) |
| `lps`, `status` | written by the check script; `status` is `ok` in all 37,620 rows |
| `FC_cplex`, `FASTCC_cplex`, `SPECTRA_cplex`, `FC_gurobi`, `FASTCC_gurobi`, `SPECTRA_gurobi` | 1 = the method kept the reaction, 0 = removed it |
| `oldFC_cplex` | the same for the earlier FasterCC release (preprint), run on 527 tasks; NaN = not run |

## How to run

```bash
python3 -I check_biomass.py            # standard library only; writes out/check_biomass_output.txt
python3 -I make_bio_wb_me.py data new_out
```

`check_biomass.py` stops with an error if its recount differs from `out/tableR23_biomass.csv`. Today it prints "re-count equals out/tableR23_biomass.csv: yes".

## Results

**Table R23 | Biomass reactions kept or removed wrongly, per method** (rows; 37,523 working, 97 blocked)

| Method | Solver | Tasks | Blocked biomass kept | Working biomass removed |
|---|---|---|---|---|
| FasterCC | CPLEX | 1,140 | 1 | 0 |
| FASTCC | CPLEX | 1,140 | 2 | 0 |
| SPECTRA | CPLEX | 1,140 | 15 | 0 |
| FasterCC (preprint) | CPLEX | 527 | 1 | 0 |
| FasterCC | Gurobi | 1,140 | 0 | 0 |
| FASTCC | Gurobi | 1,140 | 0 | 0 |
| SPECTRA | Gurobi | 1,140 | 2 | 0 |

**Why the 97 blocked biomass reactions are blocked** (from `check_biomass.py`):
- **58 are closed by their own bounds:** the task closed them (lb = ub = 0). Every method removes all 58.
- **39 have open bounds (0 to 1000)** but cannot carry flux, because reactions closed elsewhere in the task cut every route to them. All of them are `biomassreaction_N`.
  - These are the hard cases: nothing in the reaction itself says it is blocked.
  - **Every wrongly kept biomass reaction, for every method and solver, is one of these 39.**
- None of the 97 occur in tasks with 0 or 10 reactions closed: 6 occur at 50, 7 at 100, 25 at 200 and 59 at 500.

**Who keeps which:**
- **One task fools every CPLEX method:** CC_c02_l050_r09, `biomassreaction_2`. FasterCC, FASTCC, SPECTRA and the preprint FasterCC all keep it. This is the "KEPTCHK" case.
- **FASTCC with CPLEX** also keeps `biomassreaction_2` in CC_c03_l500_r07. FasterCC removes it there, and is right.
- **SPECTRA with CPLEX** keeps 15 blocked biomass reactions, in 15 tasks. With Gurobi it keeps 2, in the same two tasks as above.
- **With Gurobi**, FasterCC and FASTCC keep none.

**Message.** No method ever removes a biomass reaction that can carry flux, so none of them breaks a network's growth. The errors are all in the other direction: keeping a biomass reaction whose bounds are open but whose network can no longer feed it. FasterCC makes that error once in 2,280 runs (1,140 tasks × 2 solvers), the same once as every other CPLEX method. FASTCC makes it twice and SPECTRA 17 times.

## What this check does not show

- It covers the **cluster tasks only**. The real models are checked against FVA in `06_method_comparison_tables` (`kept_blocked_…` and `removed_consistent_…`).
- It does not say **which step of FasterCC** removed each blocked biomass reaction (dead end, exact Poolman, or LP stage). The benchmark runs did not write step records.
- To show what exact Poolman does on its own, see `toy_examples/toy_poolman.m`.
