# Unreliable heterogeneity readme

This document describes the files in the repository belonging to the article "Unreliable heterogeneity: how measurement error obscures heterogeneity in meta-analyses in psychology".

The following folders exist:
- code
- data
- figures
- manuscript
- supplement

## Code (folder)
- 01_compute_tau_fishers_z.r: code to compute tau values for Fisher's z given some tau-values of Pearson's r
- 02_functions.r: functions used for simulations
- 03_data_simulation.r: code to simulate data, parallelized.
- 04_in-texxt-values.r code to report results in the manuscript. Generates the file manuscript/in-text-values.RDS
- 05_plots.r: code to generate all plots
- 06_tables.qmd: code to create tables. Render the qmd file to create the object manuscript/tables.RDS

## Data (folder)
- means_combined.csv: simulated data. Dataframe with 15 variables and 1876 observations. Each row is the average across 10,000 replications. The first four variables are simulation results and the remainder indicate the simulation condition. The variables are:
  - intercept_b: the meta-analytic intercept
  - intercept_p: the p-value of the meta-analytic intercept
  - tau2_hat: the heterogeneity estimate from meta-analysis
  - tau2_p: the p-value of tau2_hat
  - sample_size: the sample size in the simulation condition
  - k: the number of meta-analyzed studies in the simulation condition
  - true_tau2: true heterogeneity in the simulation condition
  - mu: true average effect size in the simulation condition
  - reliability_mean: average reliability across studies in the simulation condition
  - reliability_sd: standard deviation of reliability across studies in the simulation condition
  - effect_type: Whether Pearson's correlation ("r") or Fisher's z ("r_z") in the simulation condition
  - method: whether Hedges and Vevea (HV) or Hunter and Schmidt method (HV)  in the simulation condition
  - step_length: step length of the Fisher scoring algorithm that metafor uses in estimation
  - maxiteration: max iterations fo the Fisher scoring algorithm that metafor uses in estimation
  - non_converged: number of attempted estimations that did not converge from the 10,000 repetitions.

## Figures (folder)
- manuscript (folder): figures for main text
- supplement (folder): figures for supplemental text

## Manuscript (folder)
- manuscript.qmd: quarto (dynamic) version of manuscript
- tables.RDS: RDS file with tables for main manuscript
- in-text-values.RDS: RDS file with reported results in the main manuscript
- apa.csl: apa 7th edition csl file to structure references
- unreliable-heterogeneity.bib: bib file with all references
- _extensions: template for the manuscript.qmd file to produce a nicely formatted pdf

## Supplement (folder)
- All supplements for the main text
