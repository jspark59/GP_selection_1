# GP model selection — R code for CSAM

This repository provides recovered R scripts for *Computationally Efficient Model Selection in Gaussian Process Regression*.

## Files
- `GP_Select.R` — main BF, CombinedFast, and BatchFast comparison (10 replications).
- `GP_Path.R` — exploratory beta-first versus theta-first path diagnostics.
- `GP_LOO.R` — LOO residual diagnostics and plots for Sparse Sobol10 (10 replications); reads `GP_Select.R` from the working directory.

## Usage
Keep all three R files in the same working directory. In R or RStudio:

```r
source("GP_Select.R")
```

To generate the LOO diagnostics, start a fresh R session and run:

```r
source("GP_LOO.R")
```

The main experiment defaults to n = 50, 10 replications, jitter = 1e-10, and three benchmark functions. For a quick check, change `R_pilot <- 10L` to `R_pilot <- 1L` in a working copy.

## Reproducibility notes
- These files are recovered code, **not independently revalidated** against every manuscript result. R was unavailable in the preparation environment.
- The LOO driver loads definitions from `GP_Select.R` using `eval(parse(...), envir=.GlobalEnv)` and relies on the `RESTARTABLE MONTE CARLO DRIVER` marker; do not remove that marker.
- LOO diagnostics hold covariance hyperparameters fixed at full-sample estimates; they do not repeat complete model selection.
- Before citing numerical reproducibility, execute the scripts under the manuscript settings and verify tables and figures.
