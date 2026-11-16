# FasterCC final benchmark: figure guide (16 Nov 2026)

Every figure of the reproducible pack `FasterCC_results_reproducible_16Nov`. Each entry gives the figure title and legend, a Description (what is plotted and how to read it) and a Message (what the figure shows), followed by the earlier short notes (what is plotted, message, validation). Every number below is read from the figure's own data file (`.csv` next to the figure in `out/`); nothing is estimated. Each figure now has a legend for every line, bar, point, band and reference line.

**Terms used throughout**

- *Task*: one cluster network with 0, 10, 50, 100, 200 or 500 reactions closed at random, one replicate. 19 networks × 6 levels × 10 replicates = 1,140 tasks per solver.
- *Speed-up FASTCC / FasterCC*: FASTCC's time divided by FasterCC's time in the same task and solver. Above 1, FasterCC is faster.
- *SPECTRA / FasterCC*: the same ratio for SPECTRA. Above 1, FasterCC is faster; below 1, SPECTRA is faster.
- *FVA-blocked*: a reaction that flux variability analysis says cannot carry flux. *FVA-consistent*: one that can.
- Colours are the same in every figure: **blue = FasterCC, orange = FASTCC, green = SPECTRA**. Solid line = CPLEX, dashed = Gurobi, unless the legend says otherwise.

---

## Part 1. Speed: FasterCC against FASTCC (folder 01_speed_FasterCC_vs_FASTCC)

### Figure P3 — run time by number of blocked reactions (preprint Figure 3 layout)
`figP3_runtime_by_blocked.png`

**Figure 3 | Runtime comparison of FASTCC and FasterCC across networks with increasing numbers of blocked reactions.** Median runtimes of FASTCC (orange dashed line) and FasterCC (blue line) are shown for benchmark networks containing 0–500 blocked reactions. Shaded areas indicate the interquartile range (IQR) across cluster configurations.

**Description.** Figure 3 shows how the run time of FASTCC and FasterCC changes as more reactions in the input network are blocked, with CPLEX in the left panel and Gurobi in the right. The x-axis gives the number of reactions closed at random in each cluster network (0, 10, 50, 100, 200 or 500). The y-axis gives the time per task in seconds on a logarithmic scale, so equal vertical distances mean equal ratios rather than equal seconds. Each blocked level has 190 tasks: the 19 cluster networks, from 21,200 to 183,920 reactions, each with 10 random replicates. Both methods were run on all 190 tasks with each solver, on the same input files and with no missing times, giving 1,140 tasks per solver. Each point is the median time of those 190 tasks, meaning half the tasks took less and half took more, and the line simply joins the six medians. The shaded band is the interquartile range, from the 25th to the 75th percentile. It holds the middle half of the tasks and is wide because the largest network is almost nine times the size of the smallest. FASTCC is the orange dashed line and FasterCC the blue solid line.

The shaded band shows the middle half of the 190 tasks (25th to 75th percentile). It is wide mainly because the tasks include networks of very different size, from 21,200 to 183,920 reactions, not because the times vary much between replicates of the same network.

**Message.** FASTCC’s run time grows with the number of blocked reactions, whereas FasterCC’s does not. With CPLEX, FASTCC’s median rises about fivefold, from 3,964 s with no blocked reaction to 19,630 s with 500, while FasterCC stays at 830 s and 776 s. With Gurobi, FASTCC rises from 1,429 s to 4,472 s, while FasterCC goes from 436 s to 374 s. As a result, the advantage of FasterCC grows with the number of blocked reactions. Dividing one median by the other, it goes from 4.8 to 25.3 times with CPLEX and from 3.3 to 12.0 times with Gurobi.

Computed task by task on the same network and averaged with a geometric mean, the speed-up goes from 4.0 to 25.5 times with CPLEX and from 2.7 to 11.2 times with Gurobi, and over all 1,140 tasks it is 8.11 times with CPLEX and 4.51 times with Gurobi. These per-task values are the ones to quote, because each compares the two methods on exactly the same network. This reproduces the trend of the preprint’s Figure 3 with a larger set of networks, the final FasterCC version and a second solver: FASTCC slows down as more reactions are blocked, while FasterCC stays flat.

File: `figP3_runtime_by_blocked.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** Two panels (CPLEX, Gurobi). x: reactions closed at random (0–500). y: time per task in seconds, log scale. Line with points = median over the 190 tasks of that level; shaded band = interquartile range (middle half of the tasks). Orange dashed = FASTCC, blue solid = FasterCC.
- **What “line with points = median over the 190 tasks of that level” means:**
    - *Why 190 tasks.* At each number of blocked reactions there are 19 cluster networks, each with 10 random replicates: 19 × 10 = 190 tasks. Each method was timed on all 190.
    - *The point* is the median of those 190 times: half the tasks took less, half took more. The line just joins the six points.
    - *The shaded band (IQR)* is the middle 50 % of tasks: from the 25th to the 75th percentile, which are the ranges in the tables. The band is wide because the networks range from 21,200 to 183,920 reactions.
- **How to read.** Compare how each line changes from left to right: a rising line means the method slows down as more reactions are blocked.
- **Message.** FASTCC’s time grows with the number of blocked reactions; FasterCC’s does not. CPLEX: FASTCC median 3,964 s at 0 blocked and 19,630 s at 500; FasterCC 830 s and 776 s. Gurobi: FASTCC 1,429 s → 4,472 s; FasterCC 436 s → 374 s.
- **Validation.** Same input files and solver for both methods; 1,140 tasks per solver, no missing time. Data: `figP3_runtime_by_blocked.csv`.

**Values plotted** (time per task in seconds; median, with the IQR in brackets):

| Solver | Blocked | FASTCC | FasterCC | FASTCC median ÷ FasterCC median |
|---|---|---|---|---|
| CPLEX | 0 | 3,964 (1,107–5,833) | 830 (383–1,044) | 4.8× |
| CPLEX | 10 | 4,115 (1,207–5,698) | 780 (419–979) | 5.3× |
| CPLEX | 50 | 5,546 (1,983–8,170) | 808 (385–999) | 6.9× |
| CPLEX | 100 | 6,896 (2,530–9,605) | 807 (434–1,003) | 8.5× |
| CPLEX | 200 | 10,118 (3,922–13,748) | 778 (403–1,029) | 13.0× |
| CPLEX | 500 | 19,630 (8,579–27,939) | 776 (379–960) | 25.3× |
| Gurobi | 0 | 1,429 (464–2,894) | 436 (132–908) | 3.3× |
| Gurobi | 10 | 1,532 (493–2,871) | 406 (159–850) | 3.8× |
| Gurobi | 50 | 1,841 (622–3,217) | 406 (164–944) | 4.5× |
| Gurobi | 100 | 2,054 (775–3,548) | 389 (180–840) | 5.3× |
| Gurobi | 200 | 2,541 (1,015–4,497) | 397 (152–869) | 6.4× |
| Gurobi | 500 | 4,472 (1,939–6,821) | 374 (137–819) | 12.0× |

### Figure R2 / S5 / S6 — per-task speed-up by number of blocked reactions
`figR2_S5_S6_speedup_by_blocked.png`

**Figure R2 | Per-task speed-up of FasterCC over FASTCC and SPECTRA across networks with increasing numbers of blocked reactions.** For each task, the run time of FASTCC (orange) and of SPECTRA (green) is divided by the run time of FasterCC on the same network with the same solver (left, CPLEX; right, Gurobi). Lines show the median ratio over the 190 tasks of each blocked level (0–500 blocked reactions) and shaded areas the interquartile range (IQR); values above the dashed line at 1 mean FasterCC is faster.

**Description.** Figure R2 compares the methods task by task instead of by their separate run times. For each of the 1,140 cluster tasks per solver, FASTCC’s time and SPECTRA’s time are divided by FasterCC’s time on the same network with the same solver. A value above 1 means FasterCC was faster and a value below 1 means it was slower; the grey dashed line at 1 marks equal time. The left panel is CPLEX and the right panel Gurobi. The x-axis gives the number of reactions closed at random (0, 10, 50, 100, 200 or 500). The y-axis gives the speed-up on a logarithmic scale, on which 2 and 0.5 lie at the same distance from the line at 1. Orange is FASTCC divided by FasterCC and green is SPECTRA divided by FasterCC. At each level the point is the median of the 190 task ratios (19 networks × 10 replicates), and the line joins the six medians. The shaded band shows the middle half of the 190 ratios, from the 25th to the 75th percentile. Each ratio compares two methods on the same network, so the effect of network size largely cancels out. These bands are therefore much narrower than the time bands of Figure 3, and they show how consistent the speed-up is from one network to another.

**Message.** FasterCC is faster than FASTCC in every task with CPLEX (1,140 of 1,140, by 1.5 to 59.7 times) and in 1,125 of 1,140 tasks with Gurobi, and its advantage grows with the number of blocked reactions. With CPLEX the median speed-up rises from 4.1 with no blocked reaction to 25.7 with 500, and with Gurobi from 3.1 to 11.5. The lower edge of the orange band never falls below 3.0 with CPLEX or 2.3 with Gurobi, so the gain does not depend on a few networks. The 15 Gurobi tasks where FASTCC was faster are all on the smallest network (21,200 reactions) with 0 or 10 blocked reactions. There both methods take only about 20–30 s, and the ratio is 0.75–0.99. Against SPECTRA the ratio does not change with the number of blocked reactions. With CPLEX the median stays between 0.54 and 0.58, so SPECTRA is about twice as fast as FasterCC, and FasterCC is faster in only 72 of 1,140 tasks. With Gurobi the median lies between 1.20 and 1.48, and FasterCC is faster in 822 of 1,140 tasks. The lower edge of the green Gurobi band lies just below 1 at most levels (0.90–0.98, and 1.03 with 500 blocked), because SPECTRA is still faster in 43 to 60 of the 190 tasks at each level. 238 of these 318 tasks are on the six smallest networks (up to 66,400 reactions). Blocked reactions therefore slow FASTCC down but neither FasterCC nor SPECTRA, and it is the solver that decides whether FasterCC or SPECTRA is faster.

File: `figR2_S5_S6_speedup_by_blocked.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** Two panels (CPLEX, Gurobi). x: blocked reactions. y: per-task speed-up, log scale; line = median of the 190 tasks per level, band = interquartile range. Orange = FASTCC / FasterCC, green = SPECTRA / FasterCC. Dashed grey line at 1 = equal time.
- **How to read.** Above the dashed line FasterCC is faster than the other method.
- **Message.** Against FASTCC the median speed-up rises from 4.1× (0 blocked) to 25.7× (500 blocked) with CPLEX, and from 3.1× to 11.5× with Gurobi; the whole band stays above 1. Against SPECTRA the ratio stays flat: median 0.54–0.58 with CPLEX (SPECTRA about twice as fast) and 1.2–1.5 with Gurobi (FasterCC faster).
- **Validation.** Ratios are formed task by task, so differences between networks cancel. Data: `figR2_S5_S6_speedup_by_blocked.csv`.

**Values plotted** (median ratio, with the interquartile range in brackets):

| Solver | Blocked | FASTCC ÷ FasterCC | SPECTRA ÷ FasterCC |
|---|---|---|---|
| CPLEX | 0 | 4.07 (3.05–5.66) | 0.54 (0.41–0.67) |
| CPLEX | 10 | 4.76 (3.41–6.18) | 0.56 (0.41–0.70) |
| CPLEX | 50 | 6.44 (4.52–8.27) | 0.54 (0.42–0.70) |
| CPLEX | 100 | 8.11 (6.43–10.37) | 0.55 (0.42–0.69) |
| CPLEX | 200 | 12.02 (9.90–14.87) | 0.58 (0.43–0.72) |
| CPLEX | 500 | 25.71 (20.67–31.52) | 0.55 (0.44–0.75) |
| Gurobi | 0 | 3.07 (2.33–3.48) | 1.20 (0.90–1.70) |
| Gurobi | 10 | 2.96 (2.44–3.72) | 1.30 (0.97–1.75) |
| Gurobi | 50 | 3.60 (2.90–4.39) | 1.21 (0.94–1.79) |
| Gurobi | 100 | 4.39 (3.61–5.33) | 1.31 (0.96–1.82) |
| Gurobi | 200 | 6.06 (4.63–7.67) | 1.33 (0.98–1.88) |
| Gurobi | 500 | 11.46 (8.55–15.03) | 1.48 (1.03–2.04) |

### Figure 3 — geometric-mean speed-up by number of blocked reactions
`fig3_speedup_by_level.png`

**Figure 3 (average speed-up) | Average per-task speed-up of FasterCC over FASTCC and SPECTRA by number of blocked reactions.** For each blocked level, the per-task ratios FASTCC time ÷ FasterCC time (orange) and SPECTRA time ÷ FasterCC time (green) are averaged with a geometric mean over 190 tasks (left, CPLEX; right, Gurobi); values above the dashed line at 1 mean FasterCC is faster.

**Description.** This figure shows the same ratios as Figure R2, but each point is their geometric mean instead of their median; the geometric mean is the average used in Table 1 and Table S1. For each of the 190 tasks of a blocked level (19 networks × 10 replicates), FASTCC’s time and SPECTRA’s time are divided by FasterCC’s time on the same network with the same solver. The geometric mean then averages these ratios on a logarithmic scale, so that a task in which FasterCC is twice as fast and a task in which it is twice as slow cancel out. The x-axis gives the number of reactions closed at random (0, 10, 50, 100, 200 or 500) and the y-axis the average speed-up on a logarithmic scale; the dashed line at 1 marks equal time. Orange is FASTCC divided by FasterCC and green SPECTRA divided by FasterCC; the legend in the CPLEX panel applies to both panels. No band is drawn; the spread of the ratios is shown in Figure R2.

**Message.** On average FasterCC is faster than FASTCC at every blocked level, and the gain grows with the number of blocked reactions. With CPLEX it goes from 4.0 times with no blocked reaction to 4.5, 6.2, 8.1, 12.1 and 25.5 times with 10, 50, 100, 200 and 500 blocked reactions; with Gurobi from 2.7 times to 2.9, 3.6, 4.4, 5.9 and 11.2 times. Over all 1,140 tasks the average is 8.11 times with CPLEX and 4.51 times with Gurobi. These are the speed-ups to quote, and they agree with the medians of Figure R2 (25.7 and 11.5 with 500 blocked). Against SPECTRA the average does not depend on the number of blocked reactions: it is 0.53–0.57 with CPLEX, so SPECTRA is almost twice as fast, and 1.25–1.47 with Gurobi, so FasterCC is 25–47 % faster.

File: `fig3_speedup_by_level.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** Same axes as Figure R2, but each point is the geometric mean of the 190 per-task ratios (the average used in the summary tables). Legend in the CPLEX panel applies to both panels.
- **Message.** FASTCC / FasterCC: CPLEX 4.0× at 0 blocked → 25.5× at 500; Gurobi 2.7× → 11.2×. Over all 1,140 tasks: 8.11× (CPLEX) and 4.51× (Gurobi). SPECTRA / FasterCC: CPLEX 0.53–0.57; Gurobi 1.25–1.47.
- **Validation.** Matches Table 1 and Table S1 (`table1_…`, `tableS1_…`), computed from the same `speed_tasks.tsv`.

### Figure 2 — speed-up per cluster network
`fig2_speedup_by_cluster.png`

**Figure 2 | Per-task speed-up of FasterCC over FASTCC and SPECTRA on each of the 19 cluster networks.** Box plots of the ratios FASTCC time ÷ FasterCC time (orange) and SPECTRA time ÷ FasterCC time (green) over the 60 tasks of each network (6 blocked levels × 10 replicates), with CPLEX (top) and Gurobi (bottom). Networks are ordered by size from c2 (21,200 reactions) to c20 (183,920); values above the dashed line at 1 mean FasterCC is faster.

**Description.** Figure 2 breaks the speed-up down by network instead of by number of blocked reactions. Each of the 19 cluster networks has 60 tasks: 6 blocked levels with 10 random replicates each. For every task, the time of FASTCC and the time of SPECTRA are divided by FasterCC’s time on the same task with the same solver. Each box shows the middle half of a network’s 60 ratios (25th to 75th percentile), the black bar inside it is the median, the whiskers show the spread of the remaining tasks and the dots are outlying tasks. The y-axis is logarithmic and the dashed line at 1 marks equal time. The top panel is CPLEX and the bottom panel Gurobi; orange boxes are FASTCC divided by FasterCC and green boxes SPECTRA divided by FasterCC. The networks grow in size from left to right.

**Message.** FasterCC is faster than FASTCC on every network with both solvers: the median speed-up is 4.7–10.9 times with CPLEX and 2.8–7.6 times with Gurobi. The only tasks below 1 are the 15 Gurobi tasks on the smallest network, c2. Against SPECTRA the result depends on the solver. With CPLEX SPECTRA is faster on every network, with a median ratio of 0.34–0.75. With Gurobi the median ratio is below 1 only on the five smallest networks, c2–c6 (21,200–57,360 reactions; 0.81–0.99), and above 1 on the other 14 networks (1.07–2.86, the highest on c13), so FasterCC is faster on all networks from 66,400 reactions.

File: `fig2_speedup_by_cluster.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** Box plots, one per cluster network (c2 smallest, 21,200 reactions, to c20 largest, 183,920), 60 tasks each. Top: CPLEX; bottom: Gurobi. Orange boxes = FASTCC / FasterCC, green = SPECTRA / FasterCC. Box = middle half, black bar = median, whiskers and dots = spread and outliers. Dashed line at 1 = equal time.
- **Message.** FasterCC is faster than FASTCC on every network with both solvers (median 4.7–10.9× with CPLEX, 2.8–7.6× with Gurobi). Against SPECTRA: with CPLEX SPECTRA is faster on every network (median 0.34–0.75); with Gurobi FasterCC is faster on 14 of 19 networks (medians 0.81–2.86; below 1 only on the five smallest, c2–c6, 21,200–57,360 reactions).
- **Validation.** Data: `fig2_speedup_by_cluster.csv` (median, minimum, maximum per cluster).

### Figure R1 — median time per task by network size
`figR1_median_time_by_cluster.png`

**Figure R1 | Median runtime per task of FasterCC, FASTCC and SPECTRA on each of the 19 cluster networks.** Bars give the median time over the 60 tasks of each network for FasterCC (blue), FASTCC (orange) and SPECTRA (green), with CPLEX (top) and Gurobi (bottom). Networks are ordered by size from 21,200 to 183,920 reactions; the time axis is logarithmic.

**Description.** Figure R1 shows how long each method takes on each of the 19 cluster networks. Each group of three bars is one network, labelled by its size in thousands of reactions, from 21,200 (21k) to 183,920 (183k). Each bar is the median time over that network’s 60 tasks (6 blocked levels × 10 replicates). The top panel is CPLEX and the bottom panel Gurobi; blue is FasterCC, orange FASTCC and green SPECTRA, with the legend above the panels. The y-axis is logarithmic, so equal heights mean equal ratios.

**Message.** All three methods take longer on larger networks, and FASTCC is the slowest on every network with both solvers. With CPLEX, from the smallest to the largest network, FasterCC goes from 49 s to 1,202 s, FASTCC from 250 s to 11,281 s and SPECTRA from 20 s to 691 s; SPECTRA is the fastest on every network. With Gurobi, FasterCC goes from 24 s to 1,503 s, FASTCC from 72 s to 5,001 s and SPECTRA from 19 s to 2,002 s; FasterCC and SPECTRA are close on the small networks, and FasterCC is the faster of the two on every network from 66,400 reactions.

File: `figR1_median_time_by_cluster.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** Grouped bars, one group per network (x labels give the size in thousands of reactions); top CPLEX, bottom Gurobi; y: median time per task, log scale; legend above.
- **Message.** All methods slow down with network size; FASTCC is the slowest at every size. Smallest → largest network, CPLEX: FasterCC 49 s → 1,202 s; FASTCC 250 s → 11,281 s; SPECTRA 20 s → 691 s. Gurobi: FasterCC 24 s → 1,503 s; FASTCC 72 s → 5,001 s; SPECTRA 19 s → 2,002 s.
- **Validation.** Data: `figR1_median_time_by_cluster.csv`.

### Figure S1 — median time per task by network (bars)
`figS1_time_by_cluster.png`

**Figure S1 | Median runtime per task against network size for FasterCC, FASTCC and SPECTRA.** Lines join the median times over the 60 tasks of each of the 19 cluster networks for FasterCC (blue), FASTCC (orange) and SPECTRA (green), with CPLEX (left) and Gurobi (right); the time axis is logarithmic.

**Description.** Figure S1 shows the same medians as Figure R1 as lines against network size. The x-axis gives the number of reactions in the cluster network, in thousands (21,200 to 183,920), and the y-axis the median time per task over the network’s 60 tasks, on a logarithmic scale. Each point is one network and the lines join them; blue is FasterCC, orange FASTCC and green SPECTRA. The left panel is CPLEX and the right panel Gurobi; the legend in the CPLEX panel applies to both.

**Message.** The FASTCC line lies above the two others at every size with both solvers, about ten times above FasterCC on the largest network with CPLEX (11,281 s against 1,202 s). With CPLEX the SPECTRA line lies below FasterCC at every size (691 s against 1,202 s on the largest network). With Gurobi the FasterCC and SPECTRA lines are close up to 57,360 reactions and FasterCC stays below SPECTRA from 66,400 reactions on (1,503 s against 2,002 s on the largest network). The lines show how the times grow with network size; Figure 2 gives the ratios task by task.

File: `figS1_time_by_cluster.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** Two panels (CPLEX, Gurobi). x: reactions in the cluster network (thousands). y: median time per task, log scale. One line per method (legend in the CPLEX panel).
- **Message.** Same medians as Figure R1, shown as lines; the gap between FASTCC (orange) and the two others is visible at every network size.
- **Validation.** Data: `figS1_time_by_cluster.csv`.

### Figure 4 — CPLEX against Gurobi, per method
`fig4_solver_ratio.png`

**Figure 4 | Effect of the solver on each method: runtime with CPLEX divided by runtime with Gurobi on the same task.** Lines show the geometric mean of this ratio over the 60 tasks of each cluster network for FasterCC (blue), FASTCC (orange) and SPECTRA (green); above the dashed line at 1 the method is slower with CPLEX, below it the method is faster with CPLEX.

**Description.** Figure 4 compares the two solvers for each method separately. For every task, the method’s time with CPLEX is divided by its time with Gurobi on the same input file, and the ratios of the 60 tasks of each network are averaged with a geometric mean. The x-axis gives the size of the cluster network in thousands of reactions and the y-axis the ratio on a logarithmic scale; the dashed line at 1 means the method takes the same time with both solvers. Blue is FasterCC, orange FASTCC and green SPECTRA.

**Message.** The solver does not affect the three methods in the same way, and this changes their ranking. FASTCC is 2.3–4.0 times slower with CPLEX than with Gurobi on every network. FasterCC is 0.80–2.59 times: slower with CPLEX up to 129,680 reactions, about equal at 138,720–147,760 reactions, and faster with CPLEX from 156,800 reactions (0.80–0.92). SPECTRA changes little on the small networks and is much faster with CPLEX on the large ones (0.35–0.50 from 120,640 reactions). This is why SPECTRA is faster than FasterCC with CPLEX but slower with Gurobi on the large networks (Figures 2 and 5).

File: `fig4_solver_ratio.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** x: network size. y: geometric mean of time with CPLEX / time with Gurobi for the same task, log scale. One line per method. Dashed line at 1 = same time with both solvers.
- **Message.** The solver changes the ranking. FASTCC is 2.3–4.0× slower with CPLEX on every network. FasterCC is 0.80–2.59× (slower with CPLEX up to 129,680 reactions, about equal at 138,720–147,760, faster with CPLEX from 156,800: 0.80–0.92). SPECTRA is faster with CPLEX on large networks (0.35–0.50 from cluster 13). This is why SPECTRA beats FasterCC with CPLEX but not with Gurobi.
- **Validation.** Data: `fig4_solver_ratio.csv`.

### Figure 5 / R4 — FasterCC against SPECTRA, time and relative speed
`fig5_R4_fastercc_vs_spectra.png`

**Figure 5 | Runtime comparison of FasterCC and SPECTRA across network sizes and numbers of blocked reactions, with CPLEX and Gurobi.** (A) Median runtime per task of FasterCC (blue) and SPECTRA (green) for the 19 cluster networks (21,200–183,920 reactions); solid lines CPLEX, dashed lines Gurobi. (B) Median per-task ratio of SPECTRA time to FasterCC time on the same network; above the dotted line at 1 FasterCC is faster, below it SPECTRA is faster. (C) Median runtime per task for networks with 0–500 blocked reactions. Each point is the median over the 60 tasks of one network (A, B) or the 190 tasks of one blocked level (C).

**Description.** Figure 5 compares the run time of FasterCC and SPECTRA on the same 1,140 cluster tasks per solver, run on the same input files. Panel A shows the median time per task for each of the 19 cluster networks, ordered by size from 21,200 to 183,920 reactions. Each point is the median over that network’s 60 tasks (6 blocked levels × 10 replicates). FasterCC is blue and SPECTRA green; solid lines are CPLEX and dashed lines Gurobi. The y-axis is logarithmic, so equal vertical distances mean equal ratios. Panel B divides SPECTRA’s time by FasterCC’s time task by task on the same network and shows the median of those ratios for each network (grey lines; solid CPLEX, dashed Gurobi). The dotted line at 1 marks equal time: values above it mean FasterCC was faster and values below it mean SPECTRA was faster. Panel C shows the median time per task against the number of reactions closed at random (0, 10, 50, 100, 200 or 500). Each point is the median over the 190 tasks of that level (19 networks × 10 replicates).

**Message.** Which of the two methods is faster depends on the solver, not on the number of blocked reactions. With CPLEX SPECTRA is faster on every network: the median ratio in panel B lies between 0.34 and 0.75, and on the largest network SPECTRA takes 691 s against 1,202 s for FasterCC. With Gurobi the median ratio is below 1 only on the five smallest networks (0.81–0.99) and above 1 on every network from 66,400 reactions (1.07–2.86), so FasterCC is faster on all larger networks. On the largest one FasterCC takes 1,503 s against 2,002 s for SPECTRA. The difference comes from SPECTRA, which is much faster with CPLEX than with Gurobi on large networks (691 s against 2,002 s on the largest), while FasterCC changes less between the two solvers (1,202 s against 1,503 s). Panel C shows that neither method slows down as more reactions are blocked. With CPLEX FasterCC stays between 776 and 830 s and SPECTRA between 458 and 500 s; with Gurobi FasterCC stays between 374 and 436 s and SPECTRA between 668 and 770 s. This is unlike FASTCC, whose time grows with the number of blocked reactions (Figure 3). The figure compares speed only. The two methods do not return the same networks on the large real models, where SPECTRA keeps reactions that FVA shows cannot carry flux (Figure SP4).

File: `fig5_R4_fastercc_vs_spectra.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** A: median time by network size (blue FasterCC, green SPECTRA; solid CPLEX, dashed Gurobi). B: median per-task ratio SPECTRA / FasterCC by network size (solid CPLEX, dashed Gurobi; dotted line = equal time). C: median time by number of blocked reactions (same line code as A). Every panel has its own legend.
- **Message.** With CPLEX SPECTRA is faster at every size (ratio below 1 in B). With Gurobi the median ratio is below 1 on the five smallest networks (0.81–0.99) and above 1 from cluster 7 (66,400 reactions; 1.07–2.86), so FasterCC is faster on every larger network. Neither method's time depends much on the number of blocked reactions (C).
- **Validation.** Data: `fig5_R4_fastercc_vs_spectra.csv`.

### Figure S3 — tasks where the methods keep a different number of reactions
`figS3_tasks_with_different_result.png`

**Figure S3 | Number of tasks in which FasterCC keeps a different number of reactions than FASTCC or SPECTRA, per cluster network.** Bars count, out of the 60 tasks of each network, the tasks in which the number of reactions kept by FasterCC differs from that kept by FASTCC (orange) or by SPECTRA (green), with CPLEX (left) and Gurobi (right).

**Description.** Figure S3 checks whether the three methods return consistent networks of the same size. For every task, the number of reactions kept by FasterCC is compared with the number kept by FASTCC and by SPECTRA on the same input file with the same solver. Each bar counts the tasks of one cluster network (60 tasks each) in which the numbers differ; the x-axis gives the cluster number, from 2 (smallest) to 20 (largest), and the y-axis the number of tasks. Orange bars compare FasterCC with FASTCC and green bars FasterCC with SPECTRA; the left panel is CPLEX and the right panel Gurobi.

**Message.** The three methods almost always keep the same number of reactions. FasterCC and FASTCC differ in 2 of 1,140 tasks with CPLEX (one task in cluster 2 and one in cluster 3) and in none with Gurobi. FasterCC and SPECTRA differ in 15 of 1,140 tasks with CPLEX, spread over 12 networks with at most 2 per network, and in 2 tasks with Gurobi (clusters 2 and 3). Every one of these differences is a single biomass reaction; the exact biomass check (Figure R23) shows which method is right in each case.

File: `figS3_tasks_with_different_result.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** Bars per cluster (60 tasks each): number of tasks where FasterCC keeps a different number of reactions than FASTCC (orange) or SPECTRA (green). Left CPLEX, right Gurobi.
- **Message.** FasterCC and FASTCC differ in 2 of 1,140 tasks with CPLEX (clusters 2 and 3) and in none with Gurobi. FasterCC and SPECTRA differ in 15 tasks with CPLEX and 2 with Gurobi. Each difference is one biomass reaction; which method is right is decided in Figure R23.
- **Validation.** Reaction-ID comparison: `methods_bench_16Nov_diff_ids_finished.csv` (folder 06).

### Figure SP1 — time per task, FasterCC against SPECTRA
`figSP1_task_time_fastercc_vs_spectra.png`

**Figure SP1 | Runtime of FasterCC and SPECTRA on each cluster task.** Each point is one of the 1,140 tasks per solver (left, CPLEX; right, Gurobi), with FasterCC time on the x-axis and SPECTRA time on the y-axis (logarithmic scales) and coloured by network size. Points above the solid line are tasks in which FasterCC was faster; the dotted lines mark a twofold difference.

**Description.** Figure SP1 shows every cluster task instead of medians. Each point is one task: its x position is FasterCC’s time and its y position SPECTRA’s time on the same input file with the same solver, both in seconds on logarithmic scales. The colour of the point gives the size of the network, from dark purple (21,200 reactions) to yellow (183,920), as shown in the colour bar. The solid diagonal line marks equal time; above it SPECTRA was slower and below it SPECTRA was faster. The dotted lines mark where one method takes twice as long as the other. The left panel is CPLEX and the right panel Gurobi, each with 1,140 points.

**Message.** With CPLEX almost all points lie below the solid line: SPECTRA was faster in 1,068 of 1,140 tasks and FasterCC in only 72. With Gurobi most points lie above the line: FasterCC was faster in 822 of 1,140 tasks. These are mostly the larger networks: on the 14 networks from 66,400 reactions FasterCC was faster in 732 of 840 tasks, against 90 of 300 tasks on the five smaller networks. The comparison between the two methods therefore depends on the solver and, with Gurobi, on the size of the network.

File: `figSP1_task_time_fastercc_vs_spectra.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** One point per task (1,140 per panel). x: FasterCC time; y: SPECTRA time; both log scale. Point colour = network size (colour bar). Solid line = same time; dotted lines = 2× faster or slower. Legend in each panel.
- **How to read.** Points above the solid line are tasks where SPECTRA was slower.
- **Message.** CPLEX: FasterCC faster in 72 of 1,140 tasks. Gurobi: FasterCC faster in 822 of 1,140, almost all on the larger networks (yellow/green points above the line).
- **Validation.** The 15 SPECTRA times lost to a licence failure are those of the SPEC15 reruns. Data: `figSP1_task_time_fastercc_vs_spectra.csv`.

### Figure SP2 — SPECTRA / FasterCC ratio for every network and blocked level
`figSP2_spectra_ratio_heatmap.png`

**Figure SP2 | Relative speed of SPECTRA and FasterCC for every network and number of blocked reactions.** Each cell gives the median, over 10 replicates, of SPECTRA time ÷ FasterCC time for one cluster network (rows, size in brackets) and one blocked level (columns), with CPLEX (left) and Gurobi (right). Orange cells mean FasterCC is faster, purple cells mean SPECTRA is faster, and white means equal time.

**Description.** Figure SP2 shows the SPECTRA-to-FasterCC ratio for every combination of network and number of blocked reactions. The 19 cluster networks are the rows, with their size in brackets, and the six blocked levels (0, 10, 50, 100, 200 and 500) are the columns, giving 114 cells per solver. The number in each cell is the median ratio of SPECTRA’s time to FasterCC’s time over the 10 replicates of that cell; 1.00 means equal time. The colour shows the same ratio on a log2 scale: orange when FasterCC is faster (ratio above 1), purple when SPECTRA is faster (ratio below 1) and white when both are equal. The left panel is CPLEX and the right panel Gurobi.

**Message.** With CPLEX SPECTRA is faster in 113 of the 114 cells; the only exception is cluster 11 with 500 blocked reactions (1.03), and the lowest ratio is 0.29. With Gurobi FasterCC is faster in 85 of 114 cells. The 29 cells where SPECTRA is faster are all on the six smallest networks (clusters 2–7), except one: cluster 12 with no blocked reaction (0.93). From cluster 8 (75,440 reactions) onwards every other cell is orange, up to 3.48 on cluster 13. The colour changes with the network size and the solver, not with the number of blocked reactions.

File: `figSP2_spectra_ratio_heatmap.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** Heat map: rows = 19 networks (size in brackets), columns = blocked reactions. Number in each cell = median ratio SPECTRA time / FasterCC time over its 10 replicates. Colour = log2 of that ratio: orange = FasterCC faster, purple = SPECTRA faster, white = equal.
- **Message.** CPLEX panel is purple everywhere except one cell (c11 at 500 blocked, 1.03): SPECTRA is faster. Gurobi panel: from c08 (75,440 reactions) down every cell is orange except one (c12 with 0 blocked, 0.93), up to 3.48 on c13: FasterCC is faster on large networks.
- **Validation.** Data: `figSP2_spectra_ratio_heatmap.csv`.

### Figure SP3 — real-model run times
`figSP3_real_model_times.png`

**Figure SP3 | Runtime of FasterCC, FASTCC and SPECTRA on the 58 real models.** One row per model (number of reactions in brackets), grouped by model type; FasterCC (blue circles), FASTCC (orange squares) and SPECTRA (green triangles); logarithmic time axis.

**Description.** Figure SP3 shows the run time of each method on each real model, from the small test models to the Microbiota community model (173,353 reactions) and the ME model (76,664 reactions). Each row is one model, with its number of reactions in brackets, and the rows are grouped by model type, named at the right: small test models, AGORA microbes, curated microbial models, human generic models, fetal tissues, whole-body models, Microbiota and the ME model. The x-axis gives the run time in seconds on a logarithmic scale. Blue circles are FasterCC, orange squares FASTCC and green triangles SPECTRA. Each time comes from one run per model and method. FasterCC has no point for the ME model, because its run did not finish.

**Message.** SPECTRA is the fastest method on 47 of the 57 models where all three finished, FASTCC on 9 (the small test models, two AGORA microbes, iIT341, Recon3D_mat and Recon3Model) and FasterCC on 1 (Human2). FASTCC is fast on small models but by far the slowest on large ones, and FasterCC is faster than FASTCC on 33 of the 57 models. Geometric means per group are: fetal tissues 131 s for FasterCC, 524 s for FASTCC and 13 s for SPECTRA; corrected whole-body models (aafix) 2,976 s, 117,302 s and 353 s; Microbiota 9,084 s, 342,660 s and 1,406 s. The two original whole-body models (_1_03c) take about 23 h with FasterCC, because these networks are infeasible with their own bounds and FasterCC’s precision safeguard (careful mode) then solves about 850 LPs again with single-thread primal simplex (Table R6b).

File: `figSP3_real_model_times.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** One row per real model (number of reactions in brackets), grouped by model type (label at the right). x: run time, log scale. Blue circle FasterCC, orange square FASTCC, green triangle SPECTRA; legend above.
- **Message.** SPECTRA is the fastest on almost every real model. FASTCC is fastest on small models and slowest on large ones. Geometric means per group (Table SP2): fetal tissues FasterCC 131 s, FASTCC 524 s, SPECTRA 13 s; whole-body aafix FasterCC 2,976 s, FASTCC 117,302 s, SPECTRA 353 s; Microbiota 9,084 s, 342,660 s, 1,406 s.
- **Validation and limit.** One run per model and method; times from `results_real` (Slurm elapsed times differ by a few seconds; Table R6). The two `_1_03c` whole-body models take about 23 h with FasterCC because of careful mode (Table R6b).

### Figure SP4 — real models compared with FVA
`figSP4_real_models_vs_fva.png`

**Figure SP4 | Disagreements with flux variability analysis (FVA) on the real models.** For each model where at least one method differs from FVA, bars give the number of FVA-consistent reactions the method removed (left) and the number of FVA-blocked reactions it kept (right), for FasterCC (blue), FASTCC (orange) and SPECTRA (green); symmetric logarithmic axis.

**Description.** Figure SP4 compares the output of each method with FVA of the input model, reaction by reaction counted. Only the 32 models where at least one method differs from FVA are shown, one per row. The left panel counts reactions that FVA shows can carry flux but that the method removed: this is an error, because the network loses a working reaction. The right panel counts reactions that FVA says cannot carry flux but that the method kept. The x-axis is a symmetric logarithmic scale, so that 0 is shown and large and small counts fit on one axis. Blue bars are FasterCC, orange FASTCC and green SPECTRA. FasterCC has no bar for the ME model, because its run did not finish.

**Message.** FasterCC removes a reaction that can carry flux on only 2 models, the two Harvetta models, with 3 reactions each; FASTCC removes the same 3. FASTCC removes such reactions on 10 models: Ec_iAF1260 (2), STM_v1_0 (6), iML1515, Human1 and Human2 (1 each), the four whole-body models (2–3 each) and ME (20). SPECTRA never removes one. On the other side, SPECTRA keeps FVA-blocked reactions on 7 models, all large ones (ME 2,453, Microbiota 954, Harvetta_aafix 197, Harvey_aafix 113, Harvey_1_03c 103, Harvetta_1_03c 45 and Fetus_RBC 2), and FASTCC on 6 (the four whole-body models 8–15, Microbiota 96, ME 2,453). FasterCC keeps FVA-blocked reactions on 27 models, 73 on each fetal tissue, 96 on Fetus_RBC, 57 on iIT341 and fewer on the others; these are the reactions kept on purpose by FasterCC’s guard. FasterCC is therefore the only method that never removes a working reaction except the 3 Harvetta reactions FASTCC also removes, while SPECTRA’s errors are on the side of keeping blocked reactions in the large models.

File: `figSP4_real_models_vs_fva.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** Only the models where at least one method differs from FVA. Left: FVA-consistent reactions the method removed (an error: a reaction that can carry flux is lost). Right: reactions the method kept although FVA calls them blocked. x: number of reactions, symmetric log scale (0 is shown). Bars by method; legend in each panel.
- **Message.** Left: FASTCC removes FVA-consistent reactions on Ec_iAF1260 (2), STM_v1_0 (6), iML1515, Human1, Human2 (1 each), the whole-body models and ME (20); FasterCC only on the two Harvetta models (3 each, the same 3 FASTCC removes); SPECTRA never. Right: SPECTRA keeps many FVA-blocked reactions on the large models (Microbiota 954, Harvetta_aafix 197, Harvey_aafix 113, ME 2,453); FasterCC's bars here are reactions kept by its guard, which an exact test shows are needed to keep reactions already proven consistent.
- **Validation.** FVA itself is floating-point; for ME the exact checks of 14–15 Nov are the reference. Data: `figSP4_real_models_vs_fva.csv`.

### Figure SP5 — reversible reactions left in the output
`figSP5_reversible_in_output.png`

**Figure SP5 | Share of reversible reactions in the consistent networks returned by each method.** (A) Real models: reversible reactions as a percentage of reactions kept, for SPECTRA (green triangles) and FASTCC (orange squares), which keep the input directions, and for FasterCC (blue dots), after its reorientation. (B) Cluster networks: median share over the 60 Gurobi tasks of each network, against network size.

**Description.** Figure SP5 shows how many of the reactions kept by each method are still reversible. FASTCC and SPECTRA return the reactions with their input directions, while FasterCC turns a reversible reaction into an irreversible one when it finds that the reaction can carry flux in one direction only. Panel A has one row per real model, and the x-axis gives the reversible reactions as a percentage of the reactions kept; SPECTRA (green triangles) and FASTCC (orange squares) usually overlap. Panel B shows the same percentage for the 19 cluster networks, as the median over the 60 Gurobi tasks of each network, against network size in thousands of reactions.

**Message.** FasterCC’s output has about half as many reversible reactions as FASTCC’s and SPECTRA’s. Median per model group: fetal tissues 19.4 % against 37.6 %, AGORA microbes 23.6 % against 44.3 %, curated microbial models 9.7 % against 23.2 %, human generic models 33.1 % against 46.5–46.6 %, whole-body models 26.0 % against 43.2 %, and Microbiota 24.9 % against 48.9–49.1 %. On the cluster networks FasterCC leaves 29.8–34.7 % reversible reactions against 41.7–49.5 %, and the share falls with network size for all methods. Fewer reversible reactions make later analyses cheaper, such as sampling, FVA, loop removal and MILP-based methods. Table 3 shows that the reorientations are correct: on Recon3D and Human-GEM every reaction FasterCC made irreversible can carry flux only in that direction according to FVA.

File: `figSP5_reversible_in_output.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** A: one row per real model; x = reversible reactions as % of reactions kept. Green triangle SPECTRA and orange square FASTCC (both keep the input directions), blue dot FasterCC (after its reorientation). B: the same share on the cluster networks (median of the 60 Gurobi tasks) against network size.
- **Message.** FasterCC's output has about half as many reversible reactions: fetal tissues about 19 % instead of 37 %, cluster networks 30–35 % instead of 42–50 %. Fewer reversible reactions make later sampling, FVA, loop removal and MILP-based methods cheaper.
- **Validation.** The reorientations are checked against FVA reaction by reaction in Table 3 (all 1,557 on Recon3D and 1,596 on Human-GEM confirmed). Data: `figSP5_reversible_in_output.csv`.

### Figure 4 (preprint) / R3 — FASTCORE with and without FasterCC
`figP4_R3_fastcore_preprocessing.png`

**Figure P4 | FASTCORE runtime with and without FasterCC preprocessing.** (a) Median FASTCORE time (line) and interquartile range (band) on the raw network (orange dashed) and on the network cleaned by FasterCC (blue), for core sets of 10–500 reactions; grey dotted line: cleaning plus one FASTCORE run. (b) Median per-task speed-up of FASTCORE (time without ÷ time with preprocessing) and its interquartile range; dotted line at 1: no gain.

**Description.** Figure P4 tests whether a network cleaned by FasterCC helps FASTCORE, a method that builds a context-specific model around a set of core reactions. The benchmark has 250 tasks: the cluster networks 1–10 (12,160–93,520 reactions), core sets of 10, 50, 100, 200 and 500 reactions, and five replicates each, all run with CPLEX, giving 50 tasks per core size. Panel a gives FASTCORE’s time in seconds against the core size: the line is the median over the 50 tasks and the band the middle half (25th to 75th percentile). Orange dashed is FASTCORE on the raw network and blue FASTCORE on the network cleaned by FasterCC; the grey dotted line adds the cleaning time to one FASTCORE run. Panel b divides, task by task, FASTCORE’s time on the raw network by its time on the cleaned network; above the dotted line at 1 the cleaned network is faster.

**Message.** FASTCORE is faster on the network cleaned by FasterCC in 211 of 250 tasks, and at least 10 % faster in 182. The gain grows with the size of the core set: the median per-task speed-up is 1.03 with 10 core reactions (middle half 0.90–1.27), then 1.35, 1.48, 1.48 and 1.62 with 50, 100, 200 and 500 core reactions (middle half 1.38–2.10 at 500). With 500 core reactions the median FASTCORE time falls from 254 s to 151 s. The core is kept in all 250 tasks and the F1 score does not change. Cleaning plus one FASTCORE run (480–550 s) takes longer than one FASTCORE run on the raw network (215–292 s), but the cleaning is done once per network and then serves every core set.

File: `figP4_R3_fastcore_preprocessing.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** a: FASTCORE time against core size (10–500 core reactions; 50 tasks each, clusters 1–10, CPLEX). Orange dashed = raw network, blue = network cleaned by FasterCC, line = median, band = interquartile range; grey dotted = cleaning plus one FASTCORE run. b: per-task speed-up (time without / with cleaning), median and interquartile band; dotted line = no gain.
- **Message.** FASTCORE is faster on the cleaned network in 211 of 250 tasks (at least 10 % faster in 182). The median gain grows with the core: 1.0× with 10 core reactions to 1.6× with 500. The core is kept in all 250 tasks. Cleaning costs more than one FASTCORE run (grey line), but it is done once per network and serves every core set.
- **Validation.** Data: `figP4_R3_fastcore_preprocessing.csv`; per-task table `tableR7_fastcore_by_cluster`.

### Figure R5 — correctness against independent references
`figR5_correctness.png`

**Figure R5 | Speed and correctness against independent references.** (a) Real models: SPECTRA time against FasterCC time for the 57 models with both results (logarithmic scales); point size gives the number of FVA-blocked reactions SPECTRA kept, and the dashed line marks equal time. (b) Cluster tasks with an independent FVA check: number of tasks, of 257, in which each method keeps more reactions than FVA, with CPLEX (solid bars) and Gurobi (hatched bars).

**Description.** Figure R5 puts speed and correctness side by side. In panel a each point is a real model, with FasterCC’s time on the x-axis and SPECTRA’s time on the y-axis, both in seconds on logarithmic scales. Points below the dashed line are models where SPECTRA was faster. The size of a point gives the number of reactions SPECTRA kept although FVA says they cannot carry flux (key in the legend); models with 100 or more such reactions are named. Panel b uses the 257 cluster tasks (clusters 2–6) that also have an independent FVA consistency check, run with Gurobi 12 with every LP solved to optimality. For each method it counts the tasks in which the method keeps more reactions than this reference: solid bars for CPLEX and hatched bars for Gurobi.

**Message.** In panel a most models lie below the dashed line, so SPECTRA is usually faster, and the models where SPECTRA keeps the most blocked reactions (Microbiota 954, Harvetta_aafix_E6 197, Harvey_aafix_E6 113, Harvey_1_03c 103) are among those furthest below it. For Harvey_1_03c most of the gap comes from FasterCC’s careful mode on an infeasible input (Table R6b). In panel b no method removed a reaction that the independent FVA shows can carry flux. Tasks keeping more reactions than FVA are: FasterCC 1 with CPLEX and 0 with Gurobi, FASTCC 2 and 0, SPECTRA 5 and 2. The independent check so far covers 428 of 1,420 networks and compares counts, not reaction IDs; clusters 7–20 are checked after the maintenance.

File: `figR5_correctness.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** a: real models, FasterCC time (x) against SPECTRA time (y), log scales; point size = number of FVA-blocked reactions SPECTRA kept (size key in the legend); dashed line = same time; models where SPECTRA kept 100 or more are named. b: of the 257 cluster tasks with an independent FVA check (clusters 2–6), the number where each method keeps more reactions than FVA; solid bar CPLEX, hatched bar Gurobi.
- **Message.** a: SPECTRA's speed advantage is largest exactly on the models where it keeps the most blocked reactions. b: no method removed a reaction that can carry flux; tasks keeping more than FVA: FasterCC 1 (CPLEX) / 0 (Gurobi), FASTCC 2 / 0, SPECTRA 5 / 2.
- **Limit.** The independent FVA check covers 428 of 1,420 networks so far (counts, not reaction IDs); clusters 7–20 after the maintenance.

### Figure R23 — blocked biomass reactions kept
`figR23_biomass_blocked_kept.png`

**Figure R23 | Blocked biomass reactions kept by each method in the cluster benchmark.** Bars give the number of biomass reactions that an exact LP shows cannot carry flux but that the method kept, out of 37,620 biomass rows (97 blocked), for each method and solver; light blue: the earlier FasterCC release used in the preprint (527 tasks).

**Description.** Figure R23 uses the biomass reactions as an exact test case. Every cluster task contains biomass reactions; over the 1,140 tasks there are 37,620 biomass rows, of which 37,523 can carry flux and 97 cannot, as decided by exact LPs in each task’s own input network. Each bar counts the biomass reactions that cannot carry flux but were kept by a method, for each method and solver: FasterCC (blue), FASTCC (orange) and SPECTRA (green) with CPLEX and with Gurobi. The light blue bar is the earlier FasterCC release of the preprint, run with CPLEX on 527 tasks.

**Message.** No method ever removes a biomass reaction that can carry flux: all 37,523 are kept by every method with both solvers. Biomass reactions that cannot carry flux but are kept: FasterCC 1 with CPLEX and 0 with Gurobi, FASTCC 2 and 0, and SPECTRA 15 and 2. The preprint release of FasterCC kept 1 in its 527 tasks. FasterCC’s single case (task CC_c02_l050_r09, CPLEX) is kept by every CPLEX method; it is the KEPTCHK case, fixed in FasterCC after the benchmark. FasterCC therefore has the fewest biomass errors, and SPECTRA the most.

File: `figR23_biomass_blocked_kept.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** One bar per method and solver: number of biomass reactions kept although an exact LP shows they cannot carry flux, out of 37,620 biomass rows (97 blocked). Light blue = the earlier FasterCC release (preprint version, 527 tasks).
- **Message.** No method removes a working biomass reaction. Blocked biomass kept: FasterCC 1 (CPLEX) / 0 (Gurobi), FASTCC 2 / 0, SPECTRA 15 / 2. FasterCC's one case (task CC_c02_l050_r09) is kept by every CPLEX method; it is the KEPTCHK case, fixed after the benchmark.
- **Validation.** Exact LPs in the task's own input network (`data/biocheck`).

### Figure WB — SPECTRA on the whole-body models
`figWB_spectra_wholebody.png`

**Figure WB | SPECTRA on the whole-body models compared with FVA.** For Harvetta and Harvey (aafix versions), bars give the reactions SPECTRA kept although FVA calls them blocked (green) and those it removed although FVA calls them consistent (grey), at ε = 1e-4 and ε = 1e-6, with the published SPECTRA set-up and with the coupling constraints (withC).

**Description.** Figure WB tests SPECTRA on the two corrected whole-body models, Harvetta (83,521 reactions) and Harvey (81,094), under different settings. ε is the smallest flux a reaction must reach to count as able to carry flux; the whole-body models need ε = 1e-6. “Published” is SPECTRA run as published, and “withC” is SPECTRA run with the models’ coupling constraints added. For each setting the green bar counts reactions SPECTRA kept although FVA at the same ε says they cannot carry flux, and the grey bar counts reactions it removed although FVA says they can.

**Message.** At ε = 1e-4 SPECTRA matches FVA exactly on both models. At ε = 1e-6, the value these models need, it keeps reactions that cannot carry flux: 70 on Harvetta and 62 on Harvey with the published set-up, and 430 and 558 with the coupling constraints; on Harvey with the coupling constraints it also removes 3 reactions that can carry flux. Running SPECTRA again on its own output removes further reactions (for example Harvetta at ε = 1e-6: 70,262 kept, then 70,203), so its result is not a fixed point. SPECTRA therefore matches FVA on these models only at the looser threshold.

File: `figWB_spectra_wholebody.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** For Harvetta and Harvey: green = reactions SPECTRA keeps although FVA calls them blocked; grey = reactions it removes although FVA calls them consistent. x groups: ε = 1e-4 or 1e-6; "published" = the published SPECTRA set-up, "withC" = with the coupling constraints.
- **Message.** At ε = 1e-4 SPECTRA matches FVA. At ε = 1e-6, the value the whole-body models need, it keeps 62–70 FVA-blocked reactions as published and 430–558 with the coupling constraints (Harvey withC also misses 3).
- **Validation.** Data: `tableWB_spectra_wholebody` (spectrawb_v1).

### Figure LL1 — undecided reactions per model type
`fig_ll_undec_by_type.png`

**Figure LL1 | Loop-capable reactions left undecided by the loop test, per model type.** Dots give, for each model, the percentage of loop-capable reactions the checker could not decide in the input network (grey) and in the outputs of FasterCC (blue), FASTCC (orange) and SPECTRA (green); black bars are medians; n = number of models per type.

**Description.** Figure LL1 shows how well the loop test can decide, for each method’s output, whether a reaction can carry flux outside thermodynamically infeasible loops. A loop-capable reaction is one that can carry flux in a loop; the test tries to show that it can also carry flux without a loop (loop-free) or that it cannot (loop-only), and leaves it undecided when it cannot settle the question in time. The y-axis gives the undecided reactions as a percentage of the loop-capable ones. The models are grouped by type along the x-axis, with n models per type. For each type there are four columns of dots: the input network (grey) and the outputs of FasterCC (blue), FASTCC (orange) and SPECTRA (green), each dot one model and each black bar the median.

**Message.** The loop test does not separate the three methods: within each model type their medians are almost the same. AGORA microbes have below 1 % undecided, human generic models about 11 % (10.8 % for FasterCC and SPECTRA, 11.1 % for FASTCC) and fetal tissues about 12 %. On the whole-body models about 97 % stay undecided for every method. The few points near 100 % in the fetal tissues are Fetus_RBC, where the checker reached its time limit. Loop-only reactions are therefore not where the methods differ.

File: `fig_ll_undec_by_type.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** For each model type (n = number of models), dots = one model each, coloured by which network was tested (grey = input, then FasterCC, FASTCC, SPECTRA outputs; legend below); black bar = median. y: loop-capable reactions the loop checker could not decide, %.
- **Message.** The loop test does not separate the methods: the medians are almost the same for every type (human generic about 11 %). The whole-body network stays about 97 % undecided for every method.

### Figure LL2 — how each loop-capable reaction in FasterCC's output was decided
`fig_ll_fastercc_stack.png`

**Figure LL2 | How each loop-capable reaction in FasterCC’s output was decided.** One bar per model (number of loop-capable reactions in brackets), 100 % wide: loop-free by the joint LP (step B, dark blue), loop-only in a closed island (light blue), decided by the checker (green), undecided (red).

**Description.** Figure LL2 follows each loop-capable reaction in FasterCC’s output through the steps of the loop test. Each bar is one model, with its number of loop-capable reactions in brackets, and is divided into the share of these reactions decided at each step, adding up to 100 %. Dark blue reactions were shown to carry flux without a loop by one joint LP (step B); light blue ones are loop-only because they sit in a closed island; green ones were decided one by one by the checker; red ones remain undecided.

**Message.** For most models the joint LP decides most loop-capable reactions in one step; on Recon3D, for example, it decides 6,145 of 6,968. The checker decides most of the rest, and in most models only a small share stays undecided. Two models remain almost entirely undecided: Fetus_RBC (1,044 of 1,051, where the checker reached its time limit in 1,685 of 1,722 tries) and Harvetta_aafix_E6 (24,722 of 25,484, where the checker was not run).

File: `fig_ll_fastercc_stack.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** One bar per model (number of loop-capable reactions in brackets), 100 % wide. Dark blue = shown loop-free by the joint LP (step B); light blue = loop-only, closed island; green = decided by the checker; red = undecided. Legend below.
- **Message.** For most models the joint LP decides most reactions at once (Recon3D: 6,145 of 6,968). Two models stay almost all undecided: Fetus_RBC (1,044 of 1,051; the checker reached its time limit in 1,685 of 1,722 tries) and Harvetta_aafix_E6 (24,722 of 25,484; checker not run).

### Figure LL3 — undecided reactions, FasterCC against FASTCC and SPECTRA
`fig_ll_fastercc_vs.png`

**Figure LL3 | Undecided loop-capable reactions in FasterCC’s output compared with FASTCC’s (left) and SPECTRA’s (right).** One point per model, with the number of undecided reactions plus 1 on logarithmic scales; dashed line: same number; panel titles give the number of models where FasterCC leaves fewer, the same number or more.

**Description.** Figure LL3 compares, model by model, how many loop-capable reactions the loop test leaves undecided in FasterCC’s output and in the output of the other method. Each point is one model: the x-axis gives the number undecided in FASTCC’s output (left) or SPECTRA’s output (right) and the y-axis the number in FasterCC’s output. One is added to every count so that models with none can be drawn on logarithmic scales. Points on the dashed line have the same number in both outputs; points below it are models where FasterCC leaves fewer.

**Message.** The points lie on or close to the line. Against FASTCC, FasterCC leaves fewer undecided reactions in 13 models, the same number in 26 and more in 10; against SPECTRA, fewer in 12, the same in 27 and more in 10. The loop test therefore gives almost the same result for the three methods, and loop-only reactions are not where they differ.

File: `fig_ll_fastercc_vs.png` (pdf, svg and csv in the same folder of the pack)

- **What is plotted.** One point per model; x = undecided reactions in the other method's output, y = in FasterCC's (plus 1, log scales). Dashed line = same number. Title: number of models where FasterCC has fewer, equal, more.
- **Message.** The points lie on the line: FasterCC vs FASTCC fewer 13 / equal 26 / more 10; vs SPECTRA fewer 12 / equal 27 / more 10. Loop-only reactions are not where the methods differ.

