---
output:
  pdf_document: default
  html_document: default
---

# End_of_Outbreak_Probs_Underrep

This repository contains code for reproducing the results presented in the scientific paper, "Real-time inference of the end of an outbreak: Temporally aggregated disease incidence data and under-reporting", by Isaac Ogi-Gittins, Jonathan Polonsky, Mory Keita, Steve Ahuka-Mundeke, William S Hart, Michael J Plank, Ben Lambert, Edward M Hill and Robin N Thompson.

All code written in MATLAB, compatible with version R2021b.

## Repository structure

Please find below an explainer of the directory structure within this repository.

### outputs
Contains PDFs for each figure in the manuscript.

### src

#### functions
Houses those functions called within scripts.

#### mats
MAT files with existing results. May be used to generate figure set without rerunning all simulations.

#### packages
Our analysis makes use of the [mcmcstat package](https://github.com/mjlaine/mcmcstat) for the computation of the geweke diagnostic, with associated files in the "src/pacakages/mcmcstat" directory.

#### scripts

 - **dailyAndWeeklyComps_RWD**  
Script for reproducing figures 1, 2, 3, 4, S1 and S2. 
Code for the analysis to all figures can be found in this file given above, other than figure S1, which can be found in "src/tests/convergenceOfSuperSimpleSyntheticExample.m". Note that this file may take several hours to run, if 'simSamples' variable is large (i.e. 200,000 as in the manuscript). To run the analysis for generating figures 3 and 4 (where the Gibbs-sampling method was used), you must un-comment lines 312-321.

 - **testingEmpirical_UnderrepEOO**  
For validation of the inference method outputs, a numerical simulation to estimate of end of outbreak probability in the presence of underreporting.

#### tests