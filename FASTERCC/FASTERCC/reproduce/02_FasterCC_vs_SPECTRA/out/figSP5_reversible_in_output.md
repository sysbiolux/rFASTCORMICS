**Figure SP5 | Reversible reactions left in the consistent network.** FASTCC and SPECTRA return the kept reactions with the
input bounds, so every reaction reversible in the input stays reversible. FasterCC also returns the bounds after its
reorientation (fixIrr turns backward-only reactions forward; the structure analysis makes reversible reactions irreversible
where the network and constraints allow only one direction), so fewer reversible reactions remain for downstream analysis
(sampling, FVA, loop removal, MILP methods, where each reversible reaction adds work). A: % of kept reactions that are
reversible, per real model (reversible = lb < 0 < ub; SPECTRA, open triangles, keeps the same reversible share as FASTCC and is mostly hidden under it; the same in B). B: same, median over the 60 Gurobi tasks of each cluster network
(FasterCC: res.nRev after reorientation). FasterCC's saved output model keeps the input orientation and bounds; the
reoriented counts come from its returned bounds. Data: figSP5_reversible_in_output.csv; code: make_spectra_figs.py.
