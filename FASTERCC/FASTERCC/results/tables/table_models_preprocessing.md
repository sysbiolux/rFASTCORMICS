**Table MP | Real models and FasterCC preprocessing, step by step.** One row per model: size, threshold ε, how the reactions are bounded in the input, the reactions FasterCC removes at each step (dead ends, exact Poolman test, dead ends after Poolman, LP stage), the reactions it keeps, how its preprocessing reorients reversible reactions, and any conversion of the input file before the benchmark. The 9 cluster networks without closed reactions are listed at the end.

**Description.** The table first describes each input model: its number of reactions and metabolites, the threshold ε used in the benchmark (the smallest flux a reaction must reach to count as able to carry flux), and how its reactions are bounded: reversible (lower bound below 0 and upper bound above 0), irreversible forward, irreversible backward (lower bound below 0 and upper bound at most 0) and closed (both bounds 0). The last column says whether the input file had to be converted before the methods were run, for example missing constraint senses set to equality or coupling rows moved to a separate block C. The middle columns follow the reactions through FasterCC. The preprocessing first removes dead ends: reactions of a metabolite that only one reaction uses (type I), or whose reactions are all irreversible in the same direction (type II). It then removes reactions that the exact Poolman test proves blocked, a test on the null space of the stoichiometric matrix computed in exact arithmetic with the Java engine, and the dead ends that appear after those removals. The LP stage removes the reactions that cannot reach ε in any LP, and the rest are kept as consistent. The reorientation columns count the reversible reactions that the preprocessing makes irreversible because they can carry flux in one direction only: by flipping the column when only backward flux is possible, or by changing a bound when only forward flux is possible. They also give the number of reversible reactions in the input, in the network passed to the LP stage and in FasterCC's final output. The complete file (CSV and Excel) has further columns: backward-only reactions written forward, reversible reactions flipped but still reversible, the job of each step record and FVA counts.

**Message.** Most of FasterCC's work is done before any LP is solved. Across the 54 real models with a step record (439,463 reactions in total), FasterCC removes 54,919 reactions, and 46,395 of them (84.5 %) are removed by the preprocessing: 42,909 dead ends, 3,201 reactions the exact Poolman test proves blocked and 285 dead ends that appear after those removals. The LP stage removes the other 8,524, and in 30 of the 54 models it removes none. The whole-body models follow the same pattern: on Harvey_aafix_E6 the preprocessing removes 10,272 reactions (9,316 dead ends, 939 Poolman-blocked, 17 later dead ends) and the LP stage 2,750, and on Harvetta_aafix_E6 10,620 and 2,712. Recon3D_mat, Recon3Model and the cluster networks without closed reactions are already consistent, and nothing is removed from them. The preprocessing also reduces the reversible reactions the LP stage has to handle, from 181,341 in the inputs to 96,009 in the networks passed to the LP stage. It makes 68,756 reversible reactions irreversible, 25,880 that can run only backward (column flipped) and 42,876 that can run only forward (bound changed), and writes 5,045 backward-only irreversible reactions forward. On Recon3D_mat and HumanGEM12971 these are the 1,557 and 1,596 reactions of Table 3, all confirmed by FVA.

**Validation.** The step counts are read from each run's own report files (FATE_SUMMARY, TIMELINE and REACTION_FATE, collected by `fcc_step_table.py`). For every model the steps add up to the input size, and every run's output checks are OK. The step records come from FasterCC runs of build EXA2 (jobs 16199964–16200159), not from the final benchmark runs. They keep exactly the same number of reactions as the final benchmark in all 54 models; the final version's output has slightly fewer reversible reactions on 15 models (for example Recon3D 4,253 against 4,262). No step record is available for Harvey_1_03c, Harvetta_1_03c, Microbiota and ME, so only their input description is given.


| Model | Reactions | Metabolites | ε | Reversible | Irrev. forward | Irrev. backward | Closed | Dead ends | Poolman blocked | Dead ends after Poolman | LP stage | Kept (consistent) | Reversible made irreversible | Reversible: input | at LP stage | in output | Input conversion |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **Small test** | | | | | | | | | | | | | | | | | |
| subDas | 10 | 10 | 1e-4 | 0 | 7 | 3 | 0 | 0 | 0 | 0 | 0 | 10 | 0 | 0 | 0 | 0 | missing b, csense and c set to 0, equality and 0 |
| DAS | 11 | 11 | 1e-4 | 0 | 7 | 4 | 0 | 0 | 0 | 0 | 0 | 11 | 0 | 0 | 0 | 0 | missing b, csense and c set to 0, equality and 0 |
| ecoli_core_model | 95 | 72 | 1e-4 | 46 | 49 | 0 | 0 | 0 | 0 | 0 | 8 | 87 | 20 | 46 | 26 | 26 | missing constraint senses set to equality |
| cardiac_mit_glcuptake_atpmax | 220 | 230 | 1e-4 | 98 | 110 | 3 | 9 | 117 | 27 | 0 | 0 | 76 | 22 | 98 | 24 | 24 | missing constraint senses set to equality |
| **AGORA microbes** | | | | | | | | | | | | | | | | | |
| Acidaminococcus_sp_D21 | 860 | 859 | 1e-4 | 406 | 454 | 0 | 0 | 413 | 3 | 0 | 0 | 444 | 111 | 406 | 131 | 131 | none (already in the form S·v = 0) |
| Acidaminococcus_intestini_RyC_MR95 | 1,001 | 919 | 1e-4 | 425 | 576 | 0 | 0 | 294 | 5 | 0 | 4 | 698 | 149 | 425 | 156 | 156 | none (already in the form S·v = 0) |
| Abiotrophia_defectiva_ATCC_49176 | 1,073 | 955 | 1e-4 | 475 | 598 | 0 | 0 | 303 | 3 | 0 | 0 | 767 | 196 | 475 | 169 | 169 | none (already in the form S·v = 0) |
| Acidaminococcus_fermentans_DSM_20731 | 1,102 | 1,013 | 1e-4 | 481 | 621 | 0 | 0 | 374 | 7 | 0 | 2 | 719 | 141 | 481 | 192 | 192 | none (already in the form S·v = 0) |
| Acinetobacter_junii_SH205 | 1,295 | 1,061 | 1e-4 | 494 | 801 | 0 | 0 | 383 | 20 | 10 | 0 | 882 | 195 | 494 | 207 | 207 | none (already in the form S·v = 0) |
| Acinetobacter_lwoffii_WJ10621 | 1,356 | 1,149 | 1e-4 | 537 | 819 | 0 | 0 | 431 | 7 | 0 | 0 | 918 | 191 | 537 | 224 | 224 | none (already in the form S·v = 0) |
| Acinetobacter_haemolyticus_NIPH_261 | 1,382 | 1,197 | 1e-4 | 535 | 847 | 0 | 0 | 451 | 3 | 0 | 8 | 920 | 198 | 535 | 201 | 201 | none (already in the form S·v = 0) |
| Acinetobacter_radioresistens_NIPH_2130 | 1,425 | 1,221 | 1e-4 | 596 | 829 | 0 | 0 | 454 | 4 | 0 | 0 | 967 | 203 | 596 | 246 | 246 | none (already in the form S·v = 0) |
| Acinetobacter_johnsonii_SH046 | 1,446 | 1,208 | 1e-4 | 571 | 875 | 0 | 0 | 414 | 6 | 0 | 2 | 1,024 | 204 | 571 | 247 | 247 | none (already in the form S·v = 0) |
| Acinetobacter_pittii_ANC_4052 | 1,582 | 1,300 | 1e-4 | 634 | 948 | 0 | 0 | 476 | 4 | 0 | 2 | 1,100 | 216 | 634 | 262 | 262 | none (already in the form S·v = 0) |
| Acinetobacter_calcoaceticus_PHEA_2 | 1,586 | 1,325 | 1e-4 | 643 | 943 | 0 | 0 | 462 | 5 | 0 | 0 | 1,119 | 219 | 643 | 262 | 262 | none (already in the form S·v = 0) |
| Acinetobacter_baumannii_AB0057 | 1,593 | 1,340 | 1e-4 | 660 | 933 | 0 | 0 | 484 | 5 | 0 | 0 | 1,104 | 226 | 660 | 265 | 265 | none (already in the form S·v = 0) |
| Achromobacter_xylosoxidans_A8 | 1,668 | 1,382 | 1e-4 | 638 | 1,030 | 0 | 0 | 575 | 4 | 0 | 16 | 1,073 | 227 | 638 | 250 | 250 | none (already in the form S·v = 0) |
| Achromobacter_xylosoxidans_NBRC_15126 | 1,674 | 1,369 | 1e-4 | 619 | 1,055 | 0 | 0 | 474 | 7 | 0 | 2 | 1,191 | 248 | 619 | 223 | 223 | none (already in the form S·v = 0) |
| **Curated microbial** | | | | | | | | | | | | | | | | | |
| iIT341 | 554 | 485 | 1e-4 | 183 | 368 | 3 | 0 | 49 | 4 | 0 | 67 | 434 | 113 | 183 | 55 | 50 | none (already in the form S·v = 0) |
| Sc_iND750_flux1 | 1,266 | 1,061 | 1e-4 | 436 | 830 | 0 | 0 | 379 | 36 | 74 | 146 | 631 | 162 | 436 | 106 | 97 | none (already in the form S·v = 0) |
| iMM904 | 1,577 | 1,226 | 1e-4 | 495 | 1,078 | 4 | 0 | 403 | 54 | 82 | 153 | 885 | 188 | 495 | 126 | 106 | missing constraint senses set to equality |
| Ec_iAF1260_flux1 | 2,382 | 1,668 | 1e-4 | 575 | 1,803 | 0 | 4 | 238 | 21 | 7 | 584 | 1,532 | 384 | 575 | 136 | 122 | A matrix written as S (1,668 rows); no coupling rows |
| STM_v1_0 | 2,546 | 1,802 | 1e-4 | 606 | 1,923 | 0 | 17 | 344 | 34 | 7 | 564 | 1,597 | 357 | 606 | 137 | 122 | none (already in the form S·v = 0) |
| iJO1366 | 2,583 | 1,805 | 1e-4 | 636 | 1,937 | 0 | 10 | 275 | 16 | 5 | 585 | 1,702 | 398 | 636 | 174 | 160 | missing constraint senses set to equality |
| iML1515 | 2,712 | 1,877 | 1e-4 | 663 | 2,048 | 0 | 1 | 303 | 19 | 5 | 643 | 1,742 | 397 | 663 | 183 | 169 | missing constraint senses set to equality |
| **Human generic** | | | | | | | | | | | | | | | | | |
| Recon1 | 3,742 | 2,766 | 1e-4 | 1,550 | 2,192 | 0 | 0 | 1,179 | 47 | 5 | 42 | 2,469 | 496 | 1,550 | 678 | 676 | missing constraint senses set to equality |
| Recon2 | 7,440 | 5,063 | 1e-4 | 3,270 | 3,613 | 58 | 499 | 2,053 | 60 | 3 | 7 | 5,317 | 831 | 3,270 | 1,531 | 1,531 | none (already in the form S·v = 0) |
| Recon3D_mat | 10,600 | 5,835 | 1e-4 | 5,234 | 5,331 | 35 | 0 | 0 | 0 | 0 | 0 | 10,600 | 1,557 | 5,234 | 3,677 | 3,677 | missing constraint senses set to equality |
| Recon3Model | 10,600 | 5,835 | 1e-4 | 5,234 | 5,362 | 4 | 0 | 0 | 0 | 0 | 0 | 10,600 | 1,557 | 5,234 | 3,677 | 3,677 | none (already in the form S·v = 0) |
| Human2 | 12,931 | 8,461 | 1e-4 | 5,725 | 7,198 | 0 | 8 | 1,154 | 112 | 0 | 24 | 11,641 | 2,076 | 5,725 | 3,060 | 3,058 | none (already in the form S·v = 0) |
| HumanGEM12971 | 12,971 | 8,455 | 1e-4 | 5,975 | 6,988 | 0 | 8 | 1,145 | 108 | 0 | 24 | 11,694 | 1,596 | 5,975 | 3,796 | 3,794 | none (already in the form S·v = 0) |
| HumanGEM | 13,026 | 8,365 | 1e-4 | 6,129 | 6,890 | 0 | 7 | 1,053 | 117 | 0 | 21 | 11,835 | 1,605 | 6,129 | 3,924 | 3,921 | none (already in the form S·v = 0) |
| Human1 | 13,417 | 8,461 | 1e-4 | 6,150 | 6,947 | 0 | 320 | 1,367 | 126 | 0 | 21 | 11,903 | 1,596 | 6,150 | 3,949 | 3,946 | none (already in the form S·v = 0) |
| Recon3D | 13,543 | 8,399 | 1e-4 | 6,519 | 7,024 | 0 | 0 | 1,282 | 117 | 46 | 137 | 11,961 | 1,601 | 6,519 | 4,262 | 4,253 | none (already in the form S·v = 0) |
| **Fetal tissues** | | | | | | | | | | | | | | | | | |
| Fetus_RBC | 6,519 | 4,635 | 1e-4 | 3,113 | 3,406 | 0 | 0 | 494 | 41 | 0 | 0 | 5,984 | 1,179 | 3,113 | 1,404 | 1,404 | missing constraint senses set to equality |
| Fetus_Scord | 7,786 | 5,600 | 1e-4 | 3,246 | 4,540 | 0 | 0 | 428 | 15 | 0 | 0 | 7,343 | 1,378 | 3,246 | 1,427 | 1,427 | missing constraint senses set to equality |
| Fetus_Muscle | 7,919 | 5,704 | 1e-4 | 3,269 | 4,650 | 0 | 0 | 400 | 6 | 0 | 0 | 7,513 | 1,420 | 3,269 | 1,445 | 1,445 | missing constraint senses set to equality |
| Fetus_Agland | 7,924 | 5,694 | 1e-4 | 3,249 | 4,675 | 0 | 0 | 373 | 8 | 0 | 0 | 7,543 | 1,418 | 3,249 | 1,452 | 1,452 | missing constraint senses set to equality |
| Fetus_Adipocytes | 7,972 | 5,713 | 1e-4 | 3,275 | 4,697 | 0 | 0 | 363 | 24 | 0 | 0 | 7,585 | 1,418 | 3,275 | 1,472 | 1,472 | missing constraint senses set to equality |
| Fetus_Brain | 7,986 | 5,657 | 1e-4 | 3,268 | 4,718 | 0 | 0 | 375 | 16 | 0 | 0 | 7,595 | 1,419 | 3,268 | 1,460 | 1,460 | missing constraint senses set to equality |
| Fetus_Urinarybladder | 7,990 | 5,726 | 1e-4 | 3,253 | 4,737 | 0 | 0 | 370 | 9 | 0 | 0 | 7,611 | 1,426 | 3,253 | 1,450 | 1,450 | missing constraint senses set to equality |
| Fetus_Thyroidgland | 7,992 | 5,726 | 1e-4 | 3,279 | 4,713 | 0 | 0 | 357 | 19 | 0 | 0 | 7,616 | 1,410 | 3,279 | 1,495 | 1,495 | missing constraint senses set to equality |
| Fetus_Skin | 7,999 | 5,707 | 1e-4 | 3,268 | 4,731 | 0 | 0 | 382 | 15 | 0 | 0 | 7,602 | 1,396 | 3,268 | 1,477 | 1,477 | missing constraint senses set to equality |
| Fetus_Heart | 8,012 | 5,723 | 1e-4 | 3,298 | 4,714 | 0 | 0 | 362 | 11 | 0 | 0 | 7,639 | 1,430 | 3,298 | 1,497 | 1,497 | missing constraint senses set to equality |
| Fetus_Pancreas | 8,024 | 5,736 | 1e-4 | 3,286 | 4,738 | 0 | 0 | 359 | 13 | 0 | 0 | 7,652 | 1,442 | 3,286 | 1,474 | 1,474 | missing constraint senses set to equality |
| Fetus_Colon | 8,062 | 5,760 | 1e-4 | 3,295 | 4,767 | 0 | 0 | 344 | 8 | 0 | 0 | 7,710 | 1,436 | 3,295 | 1,509 | 1,509 | missing constraint senses set to equality |
| Fetus_Lung | 8,068 | 5,763 | 1e-4 | 3,293 | 4,775 | 0 | 0 | 350 | 18 | 0 | 0 | 7,700 | 1,434 | 3,293 | 1,493 | 1,493 | missing constraint senses set to equality |
| Fetus_Spleen | 8,096 | 5,773 | 1e-4 | 3,300 | 4,796 | 0 | 0 | 356 | 9 | 0 | 0 | 7,731 | 1,458 | 3,300 | 1,479 | 1,479 | missing constraint senses set to equality |
| Fetus_Kidney | 8,136 | 5,762 | 1e-4 | 3,313 | 4,823 | 0 | 0 | 353 | 18 | 0 | 0 | 7,765 | 1,445 | 3,313 | 1,499 | 1,499 | missing constraint senses set to equality |
| Fetus_Stomach | 8,177 | 5,809 | 1e-4 | 3,347 | 4,830 | 0 | 0 | 339 | 16 | 0 | 0 | 7,822 | 1,467 | 3,347 | 1,527 | 1,527 | missing constraint senses set to equality |
| Fetus_sINE | 8,302 | 5,819 | 1e-4 | 3,362 | 4,940 | 0 | 0 | 333 | 18 | 0 | 0 | 7,951 | 1,487 | 3,362 | 1,526 | 1,526 | missing constraint senses set to equality |
| Fetus_Liver | 8,615 | 5,849 | 1e-4 | 3,584 | 5,031 | 0 | 0 | 326 | 21 | 0 | 0 | 8,268 | 1,434 | 3,584 | 1,805 | 1,805 | missing constraint senses set to equality |
| **Whole-body** | | | | | | | | | | | | | | | | | |
| Harvey_1_03c | 81,094 | – | 1e-6 | 31,713 | 45,021 | 2,898 | 1,462 | – | – | – | – | 68,233 | – | 31,713 | – | 17,682 | none (already in the form S·v = 0) |
| Harvey_aafix_E6 | 81,094 | 56,452 | 1e-6 | 31,712 | 44,945 | 2,899 | 1,538 | 9,316 | 939 | 17 | 2,750 | 68,072 | 12,487 | 31,712 | 17,930 | 17,634 | none (already in the form S·v = 0) |
| Harvetta_1_03c | 83,521 | – | 1e-6 | 33,094 | 45,975 | 2,902 | 1,550 | – | – | – | – | 70,355 | – | 33,094 | – | 18,388 | none (already in the form S·v = 0) |
| Harvetta_aafix_E6 | 83,521 | 58,851 | 1e-6 | 33,093 | 45,899 | 2,903 | 1,626 | 9,600 | 996 | 24 | 2,712 | 70,189 | 13,092 | 33,093 | 18,632 | 18,332 | none (already in the form S·v = 0) |
| **Microbiota** | | | | | | | | | | | | | | | | | |
| Microbiota | 173,353 | – | 1e-6 | 75,058 | 98,295 | 0 | 0 | – | – | – | – | 110,779 | – | 75,058 | – | 27,525 | 1 duplicate reaction ID renamed; A matrix split into S (154,801 rows) and coupling constraints C (245,811 rows) |
| **ME model** | | | | | | | | | | | | | | | | | |
| ME_dbl_E0_L1_00e_2 | 76,664 | – | 1 | 17,364 | 59,159 | 10 | 131 | – | – | – | – | run not finished | – | 17,364 | – | run not finished | 6 duplicate reaction IDs renamed; 6,086 rows with a non-zero right-hand side or an inequality moved to coupling constraints C |
| **Cluster network (no reaction closed)** | | | | | | | | | | | | | | | | | |
| Cluster01 | 12,160 | 7,395 | 1e-4 | 6,794 | 5,366 | 0 | 0 | 0 | 0 | 0 | 0 | 12,160 | 2,074 | 6,794 | 4,720 |  |  |
| Cluster02 | 21,200 | 11,670 | 1e-4 | 10,468 | 10,732 | 0 | 0 | 0 | 0 | 0 | 0 | 21,200 | 3,114 | 10,468 | 7,354 |  |  |
| Cluster03 | 30,240 | 15,945 | 1e-4 | 14,142 | 16,098 | 0 | 0 | 0 | 0 | 0 | 0 | 30,240 | 4,154 | 14,142 | 9,988 |  |  |
| Cluster04 | 39,280 | 20,220 | 1e-4 | 17,816 | 21,464 | 0 | 0 | 0 | 0 | 0 | 0 | 39,280 | 5,194 | 17,816 | 12,622 |  |  |
| Cluster05 | 48,320 | 24,495 | 1e-4 | 21,490 | 26,830 | 0 | 0 | 0 | 0 | 0 | 0 | 48,320 | 6,234 | 21,490 | 15,256 |  |  |
| Cluster06 | 57,360 | 28,770 | 1e-4 | 25,164 | 32,196 | 0 | 0 | 0 | 0 | 0 | 0 | 57,360 | 7,274 | 25,164 | 17,890 |  |  |
| Cluster07 | 66,400 | 33,045 | 1e-4 | 28,838 | 37,562 | 0 | 0 | 0 | 0 | 0 | 0 | 66,400 | 8,314 | 28,838 | 20,524 |  |  |
| Cluster08 | 75,440 | 37,320 | 1e-4 | 32,512 | 42,928 | 0 | 0 | 0 | 0 | 0 | 0 | 75,440 | 9,354 | 32,512 | 23,158 |  |  |
| Cluster09 | 84,480 | 41,595 | 1e-4 | 36,186 | 48,294 | 0 | 0 | 0 | 0 | 0 | 0 | 84,480 | 10,394 | 36,186 | 25,792 |  |  |

"–" = no FasterCC step record in the files available.
