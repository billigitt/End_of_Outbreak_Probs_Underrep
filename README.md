---
output:
  pdf_document: default
  html_document: default
---

# End_of_Outbreak_Probs_Underrep

This repository contains code for reproducing the results presented in the proposed manuscript, "Real-time inference of the end of an outbreak: Temporally aggregated disease incidence data and under-reporting".

All code written in MATLAB, compatible with version R2021b.

Code for reproducing figures 1, 2, 3, 4, S1 and S2 can be found in "src/scripts/dailyAndWeeklyComps_RWD.m".

Code for the analysis to all figures can be found in the fie given above, other than figure S1, which can be found in "src/tests/convergenceOfSuperSimpleSyntheticExample.m". Note that this file may take several hours to run, if 'simSamples' variable is large (i.e. 200,000 as in the manuscript). To run the analysis for generating figures 3 and 4 (where the Gibbs-sampling method was used), you must un-comment lines 312-321.

![fig1](outputs/Fig1.pdf)

![fig2](outputs/Fig2.pdf)

![fig3](outputs/Fig3.pdf)

![fig4](outputs/Fig4.pdf)

![figS1](outputs/FigS1.pdf)

![figS2](outputs/FigS2.pdf)
