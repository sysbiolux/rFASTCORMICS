**Table WB | SPECTRA on the corrected whole-body models against FVA.** missed = FVA-consistent but removed by SPECTRA; extra =
kept by SPECTRA but FVA-blocked (|v| >= 0.99 eps on the FVA ranges); inOK/outOK = a steady state exists in the input / in
SPECTRA's output; rerun = SPECTRA on its own output (a stable method returns the same set). 'published' = SPECTRA as
published; 'withC' = with the coupling rows. Source: spectrawb_v1_table.txt (bench2026/spectrawb_out).

| model | eps | mode | rxns | Crows | inOK | SPECTRA kept | LPs | failed LPs | sec | FVA consistent | missed | extra | outOK | rerun kept | rerun LPs | rerun failed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Harvetta_aafix | 0.0001 | published | 83,521 | 107,351 | 1 | 69,520 | 20 | 0 | 309.8 | 69,520 | 0 | 0 | 1 | 69,504 | 24 | 0 |
| Harvetta_aafix | 1e-06 | published | 83,521 | 107,351 | 1 | 70,262 | 30 | 0 | 488.0 | 70,192 | 0 | 70 | 1 | 70,203 | 30 | 0 |
| Harvetta_aafix | 1e-06 | withC | 83,521 | 107,351 | 1 | 70,622 | 58 | 1 | 958.8 | 70,192 | 0 | 430 | 1 | 70,272 | 32 | 0 |
| Harvey_aafix | 0.0001 | published | 81,094 | 103,600 | 1 | 67,436 | 20 | 0 | 215.9 | 67,436 | 0 | 0 | 1 | 67,427 | 22 | 0 |
| Harvey_aafix | 1e-06 | published | 81,094 | 103,600 | 1 | 68,134 | 20 | 1 | 247.4 | 68,072 | 0 | 62 | 1 | 68,079 | 22 | 0 |
| Harvey_aafix | 1e-06 | withC | 81,094 | 103,600 | 1 | 68,627 | 38 | 1 | 650.4 | 68,072 | 3 | 558 | 1 | 68,178 | 20 | 1 |
