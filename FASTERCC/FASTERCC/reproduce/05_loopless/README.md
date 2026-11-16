# Loop-only reactions (ledger Tables R9-R13, R18-R20; loopless figures of 15 Nov)
data/loopkit/llcheck: per model and network (INPUT, FasterCC = NEW, FASTCC, SPECTRA) the loopless checker results (cluster job outputs).
Rebuild: `python3 -I ll_permodel2.py --csv data/ll_15Nov data/loopkit/llcheck/ll_table.tsv` (table; see its header for the exact
arguments) then `python3 -I plot_loopless.py data/ll_15Nov_permodel.csv out_new`. The figures in out/ were made on 15 Nov
from data/ll_15Nov_permodel.csv (PDF and PNG).
