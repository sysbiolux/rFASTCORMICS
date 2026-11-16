# Toy examples

| File | What it does |
|---|---|
| `toy_spectra_fastcc.m` | Toy 1: builds the two-branch toy (v1 ≤ 1e-4, B leaves as 0.1 B[e]), runs FVA and SPECTRA, FASTCC and FasterCC (each again on its own output), and prints the kept sets, dead-end metabolites and whether each output is a fixed point |
| `toy_poolman.m` | Toy 2: the Poolman toy (R1–R7); exact Poolman (Java engine), FVA and the three methods |
| `check_toys_python.py` | Both toys without MATLAB: FVA (scipy HiGHS), exact Poolman (fractions) and the consistent set at ε, closed |
| `check_toys_python_output.txt` | Output of the Python check |

Run in MATLAB from this folder, with the COBRA Toolbox initialised and an LP solver set (`changeCobraSolver`):
```matlab
toy_spectra_fastcc
toy_poolman
```
FasterCC's code is not in this folder yet (see the main README). With it on the MATLAB path the scripts run it too; without it they report FasterCC as "not run". For SPECTRA, add its `src/` folder (github.com/bisect-group/SPECTRA). A method that is not on the path is skipped and reported as "not run".

The meaning of both toys and the results are explained in the main README, section "Toy examples".
