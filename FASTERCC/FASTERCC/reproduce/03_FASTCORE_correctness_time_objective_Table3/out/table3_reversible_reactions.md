**Table 3 | Number of reversible reactions by bounds, after FasterCC preprocessing, and by FVA.** Reversible by bounds: lb < 0 <
ub. Reversible by FVA: minimum flux < -1e-6 and maximum flux > 1e-6 (the same counts at 1e-9 and 1e-4). FasterCC column: the
6 Nov pipeline runs (default mode) WITHOUT FasterCC's optional reversibility check, which was not run with the current version;
every reaction its preprocessing turned irreversible is irreversible in the same direction by FVA (column 5), and the
remaining difference to FVA is what the optional check is meant to close. The preprint's FasterCC column (equal to FVA) must
be reproduced with the optional check (planned after the 16 Nov maintenance, also on iMM1415 and iMM1865). The bounds and FVA
columns reproduce the preprint exactly. Sources: data/table3/.

| model | reactions | reversible (bounds) | reversible after FasterCC preprocessing (no optional check) | preprocessing reorientations confirmed by FVA | reversible (FVA) |
|---|---|---|---|---|---|
| iMM1415 | 3,726 | 1,185 | not run with the current version |  | 363 |
| iMM1865 | 10,612 | 5,244 | not run with the current version |  | 3,531 |
| Recon3D (10,600) | 10,600 | 5,234 | 3,677 | 1557 of 1557 | 3,520 |
| Human-GEM (12,971) | 12,971 | 5,975 | 3,796 | 1596 of 1596 | 3,594 |
