**Table R22 | Objective in each method's result on the largest networks.** For each network the objective reaction is tested
in the full input and in each method's result: reactions kept; feasible with the model's bounds; largest objective flux with
the bounds opened. Job BIGOBJ 16226422 (11 Nov, all four networks with results for all methods). The whole-body _1_03c
models are infeasible with their own bounds already in the full network (BIGOBJ2 16210166: feasible only at feasibility
tolerance 1e-6), so 'max 0' there reflects the model, not the method; on Microbiota only FasterCC's network can carry
community biomass (459.2, as the full network). Source: BIGOBJ_16226422.out, BIGOBJ2_16210166.out.

| network | objective (bounds) | full input network | FasterCC | FASTCC | SPECTRA |
|---|---|---|---|---|---|
| Harvey_aafix_E6 | Whole_body_objective_rxn (1 .. 1) | 81,094; feasible yes; max 2.657 | 68,072; feasible yes; max 2.657 | 68,080; feasible yes; max 2.657 | 68,185; feasible yes; max 2.657 |
| Harvetta_aafix_E6 | Whole_body_objective_rxn (1 .. 1) | 83,521; feasible yes; max 1.275 | 70,189; feasible yes; max 1.275 | 70,197; feasible yes; max 1.275 | 70,389; feasible yes; max 1.275 |
| Microbiota | communityBiomass (0 .. 1000) | 173,353; feasible yes; max 459.2 | 110,779; feasible yes; max 459.2 | 110,846; feasible yes; max 0 | 111,704; feasible yes; max 0 |
| Harvey_1_03c | Whole_body_objective_rxn (1 .. 1) | 81,094; feasible NO; max 2.333 | 68,233; feasible NO; max 0 | 68,199; feasible NO; max 0 | 68,293; feasible NO; max 0 |
| Harvetta_1_03c | Whole_body_objective_rxn (1 .. 1) | 83,521; feasible NO; max 1.275 | 70,355; feasible NO; max 0 | 70,325; feasible NO; max 0 | 70,358; feasible NO; max 0 |
