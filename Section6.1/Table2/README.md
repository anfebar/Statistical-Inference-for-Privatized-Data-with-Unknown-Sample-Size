# Reproducing Table 2 in Section 6.1

This folder contains the MCEM materials for Table 2.

## Files

- `LRMCEM_unbounded_helper.R`: helper functions and the main MCEM implementation
- `generate_summarytable.R`: driver script to regenerate `summarytable.txt`
- `LRMCEM_analyze.R`: analysis script for the saved results
- `summarytable.txt`: final saved summary file from the original experiments

## Software

Use R with the following packages:

```r
install.packages(c("VGAM", "MASS", "dplyr", "reshape"))
```

The code also uses `parallel`, which is included with base R and does not need separate installation.

## Full run

1. Set the working directory to this folder.
2. If desired, run `generate_summarytable.R` to regenerate a new `summarytable.txt`.
3. Run `LRMCEM_analyze.R` to generate Table 2 based on `summarytable.txt`.

## Small local run

For a quick local check, edit `generate_summarytable.R` and reduce:

- the number of MCEM cycles, `num_cycles`, currently `10000`
- the number of outer iterations, `iter_vec`, currently `1:100`
- the privacy-budget values for `ep_N`, `ep_N_vec`, currently `c(0.001, 0.01, 0.1, 1, 10, Inf)`
- the summary privacy budget, `eps`, currently `1`
- the sample size, `N`, currently `1000`
- the parallel setting, `mc_cores`, currently `10`

A reasonable first local test is to reduce `num_cycles`, shorten `iter_vec`, and use a smaller `ep_N_vec`.

## Computational times

- On an Ubuntu machine with an Intel Xeon W-2145 CPU (8 physical cores, 16 logical CPUs), we estimate that one iteration run with `10000` MCEM cycles across `6` values of `ep_N` takes about `2.9` hours, assuming similar parallel resources and approximately linear scaling. The full Table 2 workflow uses `100` iterations, so users may reduce the number of iterations for smaller local checks or distribute runs across additional parallel resources for larger runs.

## HPC guidance

These analyses are highly parallelizable because the same algorithm is run repeatedly across multiple datasets, repetitions, and privacy-parameter settings. Users can therefore split jobs across those independent runs in any convenient way that matches their computing environment and available resources. In our own HPC runs, we used Slurm, but the exact parallelization strategy is left to the user.

In practice, this may require small code changes so that the jobs can be split as desired. After those jobs complete, the outputs should be collected back into the same objects or file formats expected by the downstream analysis scripts, so that the final tables, figures, and summaries can be generated in the usual way.

To generate these results, we used Purdue's Bell cluster, a community cluster that uses Slurm for job scheduling, and ran jobs in parallel using 64 CPU cores.
