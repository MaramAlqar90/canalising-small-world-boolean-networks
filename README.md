# Connectivity-dependent effects of depth-1 canalisation

Reproducibility materials for:

**Connectivity-dependent effects of depth-1 canalisation on attractor dynamics in small-world Boolean networks**

Authors: Maram Alqarni, Mark Cooper, Diane Donovan, and James Lefevre

## Contents

- `detect_attractor_canalising_reduced_threshold.cpp`  
  Compiled synchronous attractor-detection implementation.

- `simulate_canalising_SW_reduced_threshold.R`  
  Reference implementation of the directed, signed Watts--Strogatz
  simulation design with depth-1 canalising rules.

- `simulate_canalising_SW_reduced_threshold_parallel.R`  
  Driver for running the simulation workflow.

- `Analysis_and_Figures.R`  
  Fits the statistical models and generates the manuscript tables,
  figures, contrasts, cross-validation summaries, and diagnostics.

## Experimental settings

The period analysis uses:

- 8,250 original `(N, d, p, q_c, p_rep)` simulation cells;
- 100 graph--rule replicates per original cell;
- 25 initial states per graph--rule replicate;
- 825,000 graph--rule replicates in total;
- threshold perturbations redrawn for each trajectory.

The fixed-system basin-and-perturbation analysis uses:

- 7,200 fixed Boolean systems;
- 200 sampled initial states per system;
- 30 one-bit perturbation trials per system;
- a 2,000-update perturbation-return window.

## Processed analysis data

The analysis workflow uses the following processed datasets:

- `PERIOD_GRAPH_LEVEL.csv`
- `BASIN_SYSTEM_LEVEL.csv`
- `PERTURBATION_TRIAL_LEVEL.csv`
- `HAMMING_SYSTEM_TIMESTEP.csv`

These datasets contain the processed observations used to reproduce the
reported statistical analyses, tables, and figures.

## Important reproducibility note

The original production simulation code used to generate the historical
simulation outputs is no longer available. The simulation files deposited
here are a reference implementation reconstructed to match the final
Methods specification used in the manuscript and thesis.

The processed analysis datasets are the data used to reproduce the reported
statistical analyses, tables, and figures.

## Requirements

The workflow requires R and a C++ compiler compatible with Rcpp.

The analysis script automatically checks for the required R packages and
installs missing packages from CRAN.

## Run the analysis

From the repository directory, run:

```r
source("Analysis_and_Figures.R")
