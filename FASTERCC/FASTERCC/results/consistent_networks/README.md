# Consistent networks of the 58 real models

One file per model (`<model>.tsv`), one row per reaction of the **input** model, from the final benchmark (8–11 Nov 2026).
It says which reactions each method kept, so the consistent network of each method can be rebuilt from the input model.

| Column | Meaning |
|---|---|
| `rxn` | reaction ID in the input model |
| `lb_in`, `ub_in` | its bounds in the input model |
| `kept_FasterCC`, `kept_FASTCC`, `kept_SPECTRA` | 1 = the method kept the reaction (consistent), 0 = removed it, NA = the method has no finished result |
| `lb_FasterCC_out`, `ub_FasterCC_out` | the bounds in FasterCC's returned model (reversible reactions that can run one way only are written irreversible, backward-only ones forward); NA when FasterCC did not keep the reaction |

The kept counts equal the columns `kept_FasterCC`, `kept_FASTCC` and `kept_SPECTRA` of `results/tables/table_per_model.csv`
for all 58 models. The one exception is the ME model, where FasterCC's run did not finish: its column is NA, and the
table says "no result".

In the ME model, duplicate reaction IDs were renamed when the results were exported (for example `DM_leuW_tRNA_dup2`).

**Rebuilding a consistent network in MATLAB** (COBRA Toolbox), for example FasterCC's network of Recon3D:

```matlab
T = readtable('Recon3D_mat.tsv', 'FileType', 'text', 'Delimiter', '\t');
keep = T.rxn(T.kept_FasterCC == 1);
consistentModel = removeRxns(model, setdiff(model.rxns, keep));   % model = the input model
```

Source: `real_rxns.tsv`, written by `reproduce/06_method_comparison_tables/export_compare.m` from the benchmark's result
files; split into one file per model.
