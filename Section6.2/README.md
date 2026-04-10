# Reproducing Section 6.2 Results and Related Supplemental Figures

This folder contains the ATUS data application workflow used to generate the results and related supplemental figures reported for Section 6.2.

## Files

- `ATUS_code.R`: main analysis script
- `Selected_plots.R`: script to generate the reported figures
- `aux.R`: helper functions
- `slice_final.cpp`: C++ code compiled through `Rcpp`
- `data/`: input CSV files
- `figures/`: generated figure PDFs
- `output/`: location for saved `.RData` files

## Software

Use R with the following packages:

```r
install.packages(c("coda", "LaplacesDemon", "Rcpp"))
```

The script also uses the base `parallel` package that ships with R.

## Full run

1. Set the working directory to this folder.
2. Run `ATUS_code.R`.
3. This saves intermediate `.RData` files to `output/`.
4. Run `Selected_plots.R` to create the figure PDF files in `figures/`.

## Small local run

For a local check, reduce:

- the number of outer iterations, `iteration`, currently `1:10`
- the number of chains, `nch`, currently `100`
- the number of saved MCMC draws per chain, `nsave`, currently `2000`
- the retained draw indices, `draw_seq`, currently `1000 + (1:10) * 100`
- the bounded privacy-budget values, `epsilon_ss`, currently `c(1, 10)`
- the inner unbounded repetitions, `iter`, currently `1:10`
- the unbounded privacy-budget values, `epsilon_n`, currently `c(0.01, 0.1, 1, 10)`
- the parallel setting in `mclapply`, `mc.cores`, currently `15`

Because `ATUS_code.R` compiles `slice_final.cpp`, test that compilation succeeds before launching a full run.

## Computational times

- On a machine with an Intel(R) Xeon(R) W-2145 CPU @ 3.70GHz, we estimate that one outer iteration of `ATUS_code.R` takes about `19.5` hours of wall-clock time, assuming approximately linear scaling on a similar multi-core setup.
- Here, one outer iteration means that, for each value of `epsilon_S` in `{1, 10}`, the code generates one noisy bounded summary `s`, runs `100` MCMC chains of length `2000` for the bounded analysis, and runs the unbounded analysis with `100` MCMC chains of length `2000` for each of `10` realizations of `n_dp` across the `4` values of `epsilon_n`.

## HPC guidance

These analyses are highly parallelizable because the same algorithm is run repeatedly across multiple datasets, repetitions, and privacy-parameter settings. Users can therefore split jobs across those independent runs in any convenient way that matches their computing environment and available resources. In our own HPC runs, we used Slurm, but the exact parallelization strategy is left to the user.

In practice, this may require small code changes so that the jobs can be split as desired. After those jobs complete, the outputs should be collected back into the same objects or file formats expected by the downstream analysis scripts, so that the final tables, figures, and summaries can be generated in the usual way.
