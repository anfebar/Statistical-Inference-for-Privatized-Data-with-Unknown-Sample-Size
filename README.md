# Revised reproducibility materials

This folder contains a minimally reorganized version of the original reproducibility materials for the paper "Statistical Inference for Privatized Data with Unknown Sample Size." The goal is to preserve the original scripts while making the repository easier to navigate and reproduce.

## Repository layout

- `Section6.1/Table1/`: Bayesian linear-regression simulation workflow for Table 1 and related supplemental figures in Section 6.1
- `Section6.1/Table2/`: MCEM analysis for Table 2 in Section 6.1
- `Section6.1/Table3/`: parametric bootstrap analysis for Table 3 in Section 6.1
- `Section6.2/`: ATUS data application workflow for the Section 6.2 results and related supplemental figures
- `quickstart.Rmd`: reduced notebook illustrating the main workflow in a fast-running form

## Software requirements

The materials use R as the main language. Different sections require different packages.

- Table 1: `VGAM`, `MASS`, `xtable`
- Table 2: `VGAM`, `MASS`, `dplyr`, `reshape`
- Table 3: `MASS`, `parallel`, `matrixcalc`
- Section 6.2: `coda`, `LaplacesDemon`, `parallel`, `Rcpp`

For reproducibility, use a recent R installation and install the packages listed in each section-specific `README.md`.

## How to reproduce the main results

1. Choose the result of interest.
2. Change the R working directory to the corresponding subfolder.
3. Follow the instructions in that subfolder's `README.md`.
4. Use `quickstart.Rmd` first if you want a fast check of the main workflow before attempting the full analyses.

## Local versus full runs

The full analyses are computationally heavy and in some cases were originally run with parallel computing or on HPC resources. Each result-specific `README.md` includes:

- a `Full run` section describing the main scripts to run,
- a `Small local run` section pointing to the main code settings that can be reduced for a local check,
- a `Computational times` section,
- and `HPC guidance` describing how the analysis can be parallelized in a user-tailored way.
