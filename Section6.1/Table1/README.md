# Reproducing Table 1 and Related Supplemental Figures in Section 6.1

This folder contains the Bayesian linear-regression simulation workflow used for Table 1 and related supplemental figures in Section 6.1.

## Files

- `linearregression_gibbs_unbounded_fullprior.R`: main simulation script
- `UnboundedLR_Processing.R`: post-processing script for the simulation outputs
- `outputs/`: directory for saved `.RData` simulation outputs
- `figures/`: directory for PDF diagnostics written by the processing script

## Software

Use R with the following packages:

```r
install.packages(c("VGAM", "MASS", "xtable"))
```

## Full run

1. Set the working directory to this folder.
2. Run `linearregression_gibbs_unbounded_fullprior.R`.
3. The script saves generated `.RData` output files in `outputs/`.
4. Run `UnboundedLR_Processing.R` to summarize the results in `outputs/` and generate the results reported in Table 1.
5. The processing script writes figure PDF files to `figures/`.

## Small local run

For a quick local check, reduce:

- the number of outer iterations, `iteration`, currently `1:200`
- the number of MCMC cycles, `num_cycles`, currently `10000`
- the sample size, `N`, currently `1000`
- the privacy-budget values for `eps`, currently `1` for the first `100` iterations and `0.1` for the second `100`
- the privacy-budget values for `ep_N`, currently `c(10^7, 1, 0.1, 0.01, 0.001, 0.0001)`

A reasonable first local test is to reduce `num_cycles`, shorten the outer `iteration` range, and use fewer `ep_N` values. The top-level `../../quickstart.Rmd` notebook provides a faster standalone sanity check of the core privatization logic.

## Computational times

- On an Ubuntu machine with an Intel Xeon W-2145 CPU (8 physical cores, 16 logical CPUs), we estimate that one bounded run with `10000` MCMC cycles takes about `0.9` hours, while one unbounded run takes about `0.9` to `1.2` hours, depending on `ep_N`. Here, a bounded run corresponds to one combination of `eps` and `iteration`, and an unbounded run corresponds to one combination of `eps`, `iteration`, and `ep_N`. The full Table 1 simulation consists of `200` bounded runs and `1200` unbounded runs. These runs are independent and can therefore be executed in parallel. In our experiments, we used a high-performance computing environment to distribute these runs across cores, substantially reducing wall-clock time relative to serial execution.

## HPC guidance

These analyses are highly parallelizable because the same algorithm is run repeatedly across multiple datasets, repetitions, and privacy-parameter settings. Users can therefore split jobs across those independent runs in any convenient way that matches their computing environment and available resources. In our own HPC runs, we used Slurm, but the exact parallelization strategy is left to the user.

In practice, this may require small code changes so that the jobs can be split as desired. After those jobs complete, the outputs should be collected back into the same objects or file formats expected by the downstream analysis scripts, so that the final tables, figures, and summaries can be generated in the usual way.

To generate these results, we used the FSU high-performance computing (HPC) cluster, a shared research resource that uses Slurm for job scheduling, and ran jobs in parallel using 100 CPU cores.
