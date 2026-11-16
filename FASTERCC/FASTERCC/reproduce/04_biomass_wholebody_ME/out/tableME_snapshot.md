**Table ME | ME model (ME_dbl_E0_L1_00e_2, 76,658 reactions): FasterCC's state after its singleton stage against FVA and
SPECTRA.** FasterCC's run on this model had not finished (restart running); the table is its snapshot 4. FVA is a
floating-point FVA; the exact checks of 14-15 Nov (quad precision, SoPlex, Farkas) are in the report. Source: me_fva_vs_quad.tsv.

| FasterCC after the singleton stage (snapshot 4) | FVA call | SPECTRA | reactions |
|---|---|---|---|
| in A (proven consistent) | consistent | kept | 70,651 |
| still to decide (J) | blocked | kept | 2,452 |
| removed earlier | blocked | closed | 2,299 |
| still to decide (J) | blocked | closed | 1,124 |
| removed earlier | blocked | kept | 131 |
| in A (proven consistent) | blocked | kept | 1 |
