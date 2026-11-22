# CSAM Gaussian Process Model Selection — recovered R source files

These are preserved original R scripts from the earlier CSAM model-selection work. They have not been modified.

## Files
- `GP_ModelSelection_BF_CombinedFast_BatchFast_R10.R`: main 10-replication experiment for BF, COMBINED-BF (CombinedFast), and BATCH-BF (BatchFast). Base R only.
- `GP_ModelSelection_BF_CombinedFast_BatchFast_R10.R`: byte-for-byte duplicate of the `(1)` file, retained for provenance.
- `BAS_DualPath_BatchFast_Diagnostic.R`: exploratory beta-first versus theta-first path diagnostic, not the principal reported method.
- `Run_SparseSobol10_All10_LOO_Panels.R`: LOO residual diagnostics and figures for 10 sparse Sobol10 replications. Reads the `(1).R` file from the working directory.

## Suggested execution
Start R or RStudio in the directory containing these files. For the main experiment, run:

```r
source('GP_ModelSelection_BF_CombinedFast_BatchFast_R10.R')
```

To generate LOO panels, run in a fresh R session:

```r
source('Run_SparseSobol10_All10_LOO_Panels.R')
```

The main experiment defaults to 10 replicates at n=50, jitter=1e-10, three benchmark functions. The scripts create output directories and can be computationally intensive. For a smoke test, change `R_pilot <- 10L` to `R_pilot <- 1L` in a *working copy* before running.

## Important repository notes
- These are recovered **original** files, not an independently revalidated release. The current environment does not have R installed, so they have not been executed or checked against every table in the submitted manuscript.
- The LOO driver uses `eval(parse(...), envir=.GlobalEnv)` to load definitions from the main script and expects its `RESTARTABLE MONTE CARLO DRIVER` marker to remain unchanged.
- For a clean public GitHub repository, consider retaining one canonical copy of the duplicated main script and changing the LOO script reference accordingly, after testing.
- Fixed covariance hyperparameters are used during LOO calculation: this is not a full leave-one-out refitting of model selection.
- Before public release, verify the numerical results and plots using the same R version and computational settings used for the manuscript.
