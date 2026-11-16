**Table 3 check | FasterCC preprocessing against FVA, reaction by reaction.** Every reaction FasterCC reoriented is found by ID
in the FVA file of the same model (same reaction count), and every reversible reaction it made irreversible can carry flux
only in that direction by FVA. Source: data/table3/.

| model | reactions reoriented by FasterCC (any kind) | IDs not found in the FVA file | reversible -> irreversible | same direction in FVA | other direction in FVA | left reversible by preprocessing minus reversible in FVA |
|---|---|---|---|---|---|---|
| Recon3D (10,600) | 1,989 | 0 | 1,557 | 1,557 | 0 | 157 |
| Human-GEM (12,971) | 2,026 | 0 | 1,596 | 1,596 | 0 | 202 |
