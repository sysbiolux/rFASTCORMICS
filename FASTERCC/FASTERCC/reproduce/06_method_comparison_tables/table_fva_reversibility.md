**Table 3 (extended) | Reversible reactions by FVA against FasterCC, and the gap left for the optional reversibility check.** For each model: the FVA class of every reaction (reversible, forward only, backward only, blocked), the reversible reactions by the input bounds, the reversible reactions FasterCC's preprocessing made irreversible and whether FVA agrees with each direction, the reversible reactions left in FasterCC's consistent network, and the gap to FVA. FasterCC's optional reversibility check was not run with the final version.

**Description.** FVA gives, for every reaction, its smallest and largest possible flux. A reaction is reversible by FVA if its minimum is below −1e-6 and its maximum above 1e-6, forward only or backward only if it can carry flux in one direction only, and blocked if it can carry none (the same threshold as Table 3). The FasterCC columns come from the step-record runs in default mode, the runs used for Table 3. "Made irreversible by preprocessing" counts the reversible reactions that FasterCC's preprocessing turned into irreversible ones. "Same direction in FVA" checks each of them by reaction ID: a reaction FasterCC turned forward must have no backward flux in FVA, and one it flipped must have no forward flux. "Reversible in FasterCC's network" counts the kept reactions that are still reversible after the preprocessing. The gap is that count minus the reversible reactions by FVA: the reactions FasterCC leaves reversible although FVA says they can run in one direction only. Closing this gap is the job of FasterCC's optional reversibility check (5th argument revopt = 1), which was not run with the final version because it adds run time and is not needed for consistency. The last columns give the reversible reactions in the final benchmark's output and the FVA file used.

**Message (models with an FVA file so far).** Every reaction FasterCC's preprocessing made irreversible runs in that direction only by FVA: 1,557 of 1,557 on Recon3D_mat and 1,596 of 1,596 on HumanGEM12971, with none in the other direction. The preprocessing therefore never removes a direction that FVA allows. It does not find all one-directional reactions: FasterCC leaves 3,677 reversible reactions on Recon3D_mat against 3,520 by FVA (gap 157), and 3,796 on HumanGEM12971 against 3,594 (gap 202). These 157 and 202 reactions are what the optional reversibility check would close. For iMM1415 and iMM1865, the two other models of the preprint's Table 3, only FVA is available (363 and 3,531 reversible); FasterCC has not been run on them with the final version.

**Status.** The other 63 models and networks have FasterCC step records but no FVA file in the files available here, so their FVA columns and gap are empty. The FVA results exist on the HPC; once they are copied, the same script (`make_table_fva_reversibility.py`) fills every row.


| Model | Reactions | FVA: reversible | forward only | backward only | blocked | Reversible by bounds | Made irreversible by FasterCC preprocessing | same direction in FVA | other direction in FVA | Reversible in FasterCC network | Gap to FVA (optional check) | Reversible in final output |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **Small test** | | | | | | | | | | | | |
| subDas | 10 | – | – | – | – | 0 | 0 | – | – | 0 | – | 0 |
| DAS | 11 | – | – | – | – | 0 | 0 | – | – | 0 | – | 0 |
| ecoli_core_model | 95 | – | – | – | – | 46 | 20 | – | – | 26 | – | 26 |
| cardiac_mit_glcuptake_atpmax | 220 | – | – | – | – | 98 | 22 | – | – | 24 | – | 24 |
| **AGORA microbes** | | | | | | | | | | | | |
| Acidaminococcus_sp_D21 | 860 | – | – | – | – | 406 | 111 | – | – | 131 | – | 131 |
| Acidaminococcus_intestini_RyC_MR95 | 1,001 | – | – | – | – | 425 | 149 | – | – | 156 | – | 156 |
| Abiotrophia_defectiva_ATCC_49176 | 1,073 | – | – | – | – | 475 | 196 | – | – | 169 | – | 169 |
| Acidaminococcus_fermentans_DSM_20731 | 1,102 | – | – | – | – | 481 | 141 | – | – | 192 | – | 192 |
| Acinetobacter_junii_SH205 | 1,295 | – | – | – | – | 494 | 195 | – | – | 207 | – | 207 |
| Acinetobacter_lwoffii_WJ10621 | 1,356 | – | – | – | – | 537 | 191 | – | – | 224 | – | 224 |
| Acinetobacter_haemolyticus_NIPH_261 | 1,382 | – | – | – | – | 535 | 198 | – | – | 201 | – | 201 |
| Acinetobacter_radioresistens_NIPH_2130 | 1,425 | – | – | – | – | 596 | 203 | – | – | 246 | – | 246 |
| Acinetobacter_johnsonii_SH046 | 1,446 | – | – | – | – | 571 | 204 | – | – | 247 | – | 247 |
| Acinetobacter_pittii_ANC_4052 | 1,582 | – | – | – | – | 634 | 216 | – | – | 262 | – | 262 |
| Acinetobacter_calcoaceticus_PHEA_2 | 1,586 | – | – | – | – | 643 | 219 | – | – | 262 | – | 262 |
| Acinetobacter_baumannii_AB0057 | 1,593 | – | – | – | – | 660 | 226 | – | – | 265 | – | 265 |
| Achromobacter_xylosoxidans_A8 | 1,668 | – | – | – | – | 638 | 227 | – | – | 250 | – | 250 |
| Achromobacter_xylosoxidans_NBRC_15126 | 1,674 | – | – | – | – | 619 | 248 | – | – | 223 | – | 223 |
| **Curated microbial** | | | | | | | | | | | | |
| iIT341 | 554 | – | – | – | – | 183 | 113 | – | – | 55 | – | 50 |
| Sc_iND750_flux1 | 1,266 | – | – | – | – | 436 | 162 | – | – | 105 | – | 97 |
| iMM904 | 1,577 | – | – | – | – | 495 | 188 | – | – | 121 | – | 106 |
| Ec_iAF1260_flux1 | 2,382 | – | – | – | – | 575 | 385 | – | – | 134 | – | 122 |
| STM_v1_0 | 2,546 | – | – | – | – | 606 | 358 | – | – | 135 | – | 122 |
| iJO1366 | 2,583 | – | – | – | – | 636 | 399 | – | – | 172 | – | 160 |
| iML1515 | 2,712 | – | – | – | – | 663 | 398 | – | – | 178 | – | 169 |
| **Human generic** | | | | | | | | | | | | |
| Recon1 | 3,742 | – | – | – | – | 1,550 | 496 | – | – | 677 | – | 676 |
| Recon2 | 7,440 | – | – | – | – | 3,270 | 831 | – | – | 1,531 | – | 1,531 |
| Recon3D_mat | 10,600 | 3,520 | 6,595 | 485 | 0 | 5,234 | 1,557 | 1,557 | 0 | 3,677 | 157 | 3,677 |
| Recon3Model | 10,600 | – | – | – | – | 5,234 | 1,557 | – | – | 3,677 | – | 3,677 |
| Human2 | 12,931 | – | – | – | – | 5,725 | 2,076 | – | – | 3,060 | – | 3,058 |
| HumanGEM12971 | 12,971 | 3,594 | 7,481 | 619 | 1,277 | 5,975 | 1,596 | 1,596 | 0 | 3,796 | 202 | 3,794 |
| HumanGEM | 13,026 | – | – | – | – | 6,129 | 1,605 | – | – | 3,924 | – | 3,921 |
| Human1 | 13,417 | – | – | – | – | 6,150 | 1,596 | – | – | 3,949 | – | 3,946 |
| Recon3D | 13,543 | – | – | – | – | 6,519 | 1,601 | – | – | 4,262 | – | 4,253 |
| **Fetal tissues** | | | | | | | | | | | | |
| Fetus_RBC | 6,519 | – | – | – | – | 3,113 | 1,179 | – | – | 1,404 | – | 1,404 |
| Fetus_Scord | 7,786 | – | – | – | – | 3,246 | 1,378 | – | – | 1,427 | – | 1,427 |
| Fetus_Muscle | 7,919 | – | – | – | – | 3,269 | 1,420 | – | – | 1,445 | – | 1,445 |
| Fetus_Agland | 7,924 | – | – | – | – | 3,249 | 1,418 | – | – | 1,452 | – | 1,452 |
| Fetus_Adipocytes | 7,972 | – | – | – | – | 3,275 | 1,418 | – | – | 1,472 | – | 1,472 |
| Fetus_Brain | 7,986 | – | – | – | – | 3,268 | 1,419 | – | – | 1,460 | – | 1,460 |
| Fetus_Urinarybladder | 7,990 | – | – | – | – | 3,253 | 1,426 | – | – | 1,450 | – | 1,450 |
| Fetus_Thyroidgland | 7,992 | – | – | – | – | 3,279 | 1,410 | – | – | 1,495 | – | 1,495 |
| Fetus_Skin | 7,999 | – | – | – | – | 3,268 | 1,396 | – | – | 1,477 | – | 1,477 |
| Fetus_Heart | 8,012 | – | – | – | – | 3,298 | 1,430 | – | – | 1,497 | – | 1,497 |
| Fetus_Pancreas | 8,024 | – | – | – | – | 3,286 | 1,442 | – | – | 1,474 | – | 1,474 |
| Fetus_Colon | 8,062 | – | – | – | – | 3,295 | 1,436 | – | – | 1,509 | – | 1,509 |
| Fetus_Lung | 8,068 | – | – | – | – | 3,293 | 1,434 | – | – | 1,493 | – | 1,493 |
| Fetus_Spleen | 8,096 | – | – | – | – | 3,300 | 1,458 | – | – | 1,479 | – | 1,479 |
| Fetus_Kidney | 8,136 | – | – | – | – | 3,313 | 1,445 | – | – | 1,499 | – | 1,499 |
| Fetus_Stomach | 8,177 | – | – | – | – | 3,347 | 1,467 | – | – | 1,527 | – | 1,527 |
| Fetus_sINE | 8,302 | – | – | – | – | 3,362 | 1,487 | – | – | 1,526 | – | 1,526 |
| Fetus_Liver | 8,615 | – | – | – | – | 3,584 | 1,434 | – | – | 1,805 | – | 1,805 |
| **Whole-body** | | | | | | | | | | | | |
| Harvey_1_03c | 81,094 | – | – | – | – | 31,713 | no step record | – | – | no step record | – | 17,682 |
| Harvey_aafix_E6 | 81,094 | – | – | – | – | 31,712 | no step record | – | – | no step record | – | 17,634 |
| Harvetta_1_03c | 83,521 | – | – | – | – | 33,094 | no step record | – | – | no step record | – | 18,388 |
| Harvetta_aafix_E6 | 83,521 | – | – | – | – | 33,093 | no step record | – | – | no step record | – | 18,332 |
| **Microbiota** | | | | | | | | | | | | |
| Microbiota | 173,353 | – | – | – | – | 75,058 | no step record | – | – | no step record | – | 27,525 |
| **ME model** | | | | | | | | | | | | |
| ME_dbl_E0_L1_00e_2 | 76,664 | – | – | – | – | 17,364 | no step record | – | – | no step record | – | run not finished |
| **Cluster network (no reaction closed)** | | | | | | | | | | | | |
| Cluster01 | 12,160 | – | – | – | – | 6,794 | 2,074 | – | – | 4,720 | – | – |
| Cluster02 | 21,200 | – | – | – | – | 10,468 | 3,114 | – | – | 7,354 | – | – |
| Cluster03 | 30,240 | – | – | – | – | 14,142 | 4,154 | – | – | 9,988 | – | – |
| Cluster04 | 39,280 | – | – | – | – | 17,816 | 5,194 | – | – | 12,622 | – | – |
| Cluster05 | 48,320 | – | – | – | – | 21,490 | 6,234 | – | – | 15,256 | – | – |
| Cluster06 | 57,360 | – | – | – | – | 25,164 | 7,274 | – | – | 17,890 | – | – |
| Cluster07 | 66,400 | – | – | – | – | 28,838 | 8,314 | – | – | 20,524 | – | – |
| Cluster08 | 75,440 | – | – | – | – | 32,512 | 9,354 | – | – | 23,158 | – | – |
| Cluster09 | 84,480 | – | – | – | – | 36,186 | 10,394 | – | – | 25,792 | – | – |
| **Preprint Table 3 (not in the benchmark)** | | | | | | | | | | | | |
| iMM1415 | 3,726 | 363 | 1,825 | 244 | 1,294 | 1,185 | not run with the final version | – | – | not run with the final version | – | – |
| iMM1865 | 10,612 | 3,531 | 6,595 | 486 | 0 | 5,244 | not run with the final version | – | – | not run with the final version | – | – |

"–" = FVA file not available here, or not applicable.
