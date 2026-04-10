# Reproducing Table 3 in Section 6.1

This folder contains the parametric bootstrap materials used for Table 3.

## Files

- `linearregression_parametric_bootstrap.R`: main simulation script
- `DPBootsLMFunctions.R`: bootstrap helper functions
- `Bootstrap_Processing.R`: post-processing script

## Software

Use R with the following packages:

```r
install.packages(c("VGAM", "MASS", "parallel", "matrixcalc", "xtable"))
```

## Full run

1. Set the working directory to this folder.
2. Run `linearregression_parametric_bootstrap.R`.
3. This writes `bootstrap_outputs.csv`.
4. Run `Bootstrap_Processing.R` to summarize the results and reproduce Table 3.

## Small local run

For a local check, reduce:

- the number of outer iterations in the outer `mclapply`, currently `1:100`
- the privacy-parameter matrix, `Eps`, currently containing the full set of `(ep_N, eps)` combinations:
  - bounded-like row `(1e7, 1e7)`
  - rows with `eps = 1` and `ep_N in c(1e7, 1, 0.1, 0.01, 0.001, 0.0001)`
  - rows with `eps = 0.1` and `ep_N in c(1e7, 1, 0.1, 0.01, 0.001, 0.0001)`
- the sample size, `N`, currently `1000`
- the inner parallel setting in `mclapply`, `mc.cores`, currently `3`
- the outer parallel setting in `mclapply`, `mc.cores`, currently `5`

Running a single outer iteration and a reduced `Eps` matrix is a reasonable first test.

## Computational times

- In an Intel(R) Xeon(R) W-2145 CPU @ 3.70GHz processor with 16 cores available for parallel computing, running `linearregression_parametric_bootstrap.R` with `100` outer iterations, `10000` noise draws, and `10000` bootstrap draws completed in about `13.8` minutes of wall-clock time (`elapsed = 829.331` seconds in `linearregression_parametric_bootstrap.Rout`). This run used nested parallelism on the 16-core machine described above.

## HPC guidance

These analyses are highly parallelizable because the same algorithm is run repeatedly across multiple datasets, repetitions, and privacy-parameter settings. Users can therefore split jobs across those independent runs in any convenient way that matches their computing environment and available resources. In our own HPC runs, we used Slurm, but the exact parallelization strategy is left to the user.

In practice, this may require small code changes so that the jobs can be split as desired. After those jobs complete, the outputs should be collected back into the same objects or file formats expected by the downstream analysis scripts, so that the final tables, figures, and summaries can be generated in the usual way.
