# Sparse Sobol10: all 10 replications x (BF, CombinedFast, BatchFast)
# Figure for the paper: 3 representative replications x 3 methods.
# Appendix: remaining 7 replications x 3 methods (large-format and 4+3 print pages).
# Requires GP_ModelSelection_BF_CombinedFast_BatchFast_R10.R in the working directory.
# Run in RStudio: source('Run_SparseSobol10_All10_LOO_Panels.R')
# IMPORTANT: LOO uses fitted covariance hyperparameters and selected variables
# from the full data, i.e. it does not re-run model selection per left-out case.

original <- 'GP_ModelSelection_BF_CombinedFast_BatchFast_R10.R'
if (!file.exists(original)) stop('Missing original R script: ', original)
lines <- readLines(original, warn = FALSE)
cut <- grep('RESTARTABLE MONTE CARLO DRIVER', lines, fixed = TRUE)[1]
if (is.na(cut)) stop('Could not locate driver boundary in original script.')
# Read only definitions (the original 10x3 Monte Carlo driver is NOT executed).
eval(parse(text = paste(lines[seq_len(cut - 1L)], collapse='\n')), envir=.GlobalEnv)

outdir_all <- 'SparseSobol10_LOO_10rep_results'
dir.create(outdir_all, showWarnings=FALSE, recursive=TRUE)
algorithms_all <- c('BF','COMBINED-BF','BATCH-BF')
labels_all <- c('BF','CombinedFast','BatchFast')
colors_all <- c('#303030','#C5574A','#376EAA')
pch_all <- c(16L,17L,15L)

calculate_loo <- function(fit, dat, method_name, rep_id) {
  X0 <- as.matrix(dat$X); y0 <- as.numeric(dat$y)
  F0 <- as.matrix(make_F(dat$X, fit$b))
  R0 <- cov_matrix(dat$X, dat$X, fit$t, fit$ell)
  diag(R0) <- diag(R0) + main_jitter
  U0 <- chol(R0)
  W0 <- chol2inv(U0)
  WF <- W0 %*% F0
  A0 <- crossprod(F0, WF)
  P0 <- W0 - WF %*% solve(A0, t(WF))
  P0 <- (P0 + t(P0))/2
  if (any(!is.finite(diag(P0))) || any(diag(P0) <= 0))
    stop('Invalid LOO residual precision diagonal')
  e <- as.vector(P0 %*% y0)/diag(P0)
  loo_fitted <- y0 - e
  train_fitted <- as.vector(predict_gp(fit, dat$X))
  rmspe <- sqrt(mean((predict_gp(fit, dat$Xtest)-dat$ytest)^2))
  list(data=data.frame(replication=rep_id,algorithm=method_name,
                       index=seq_along(y0),fitted=train_fitted,
                       loo_fitted=loo_fitted,loo_error=e,
                       stringsAsFactors=FALSE),
       stats=data.frame(replication=rep_id,algorithm=method_name,
                        Sbeta=paste(fit$b,collapse=','),
                        Stheta=paste(fit$t,collapse=','),
                        BIC=fit$BIC,RMSPE=rmspe,
                        LOO_RMSE=sqrt(mean(e^2)),
                        MaxAbs_LOO=max(abs(e)),
                        stringsAsFactors=FALSE))
}

run_one_rep_alg <- function(rep_id, alg, dat) {
  cfg$n <<- 50L; cfg$d <<- 10L; cfg$jitter <<- main_jitter
  cfg$selection_mode <<- if (startsWith(alg,'BATCH-')) 'batch' else
    if (startsWith(alg,'COMBINED-')) 'combined' else 'reference'
  X <<- dat$X; y <<- dat$y; Xtest <<- dat$Xtest; ytest <<- dat$ytest
  allvars <<- paste0('x',1:10)
  cache <<- new.env(parent=emptyenv())
  fit_calls <<- 0L; trace_rows <<- list()
  initial <- fit_model(character(0),allvars)
  if (is.null(initial)) stop('Initial model fit failed')
  result <- run_search('BF',initial)
  calculate_loo(result$fit,dat,alg,rep_id)
}

# Checkpoint one method at a time; rerunning resumes already-computed results.
# No fake or approximate residual data are used.
for (rep_id in c(4L, setdiff(1:10,4L))) {
  dat <- make_sample('sparse_sobol10',50L,rep_id)
  for (alg in algorithms_all) {
    file <- file.path(outdir_all,sprintf('rep%02d_%s.rds',rep_id,alg))
    if (file.exists(file)) {
      x <- readRDS(file)
      if (!identical(x$stats$replication,as.integer(rep_id)) ||
          !identical(x$stats$algorithm,alg))
        stop('Checkpoint content mismatch: ',file)
    } else {
      cat(sprintf('Computing replication %02d / %s ...\n',rep_id,alg))
      x <- run_one_rep_alg(rep_id,alg,dat)
      saveRDS(x,file)
    }
    cat(sprintf('rep=%02d %-12s BIC=%10.5f RMSPE=%.5f LOO_RMSE=%.5f\n',
                rep_id,alg,x$stats$BIC,x$stats$RMSPE,x$stats$LOO_RMSE))
  }
}

all_results <- lapply(1:10,function(rep_id) {
  lapply(algorithms_all,function(alg) {
    readRDS(file.path(outdir_all,sprintf('rep%02d_%s.rds',rep_id,alg)))
  })
})
all_flat <- do.call(c, all_results)
points_df <- do.call(rbind,lapply(all_flat,`[[`,'data'))
stats_df <- do.call(rbind,lapply(all_flat,`[[`,'stats'))
rownames(points_df) <- NULL; rownames(stats_df) <- NULL
write.csv(points_df,file.path(outdir_all,'LOO_all_30_models.csv'),row.names=FALSE)
write.csv(stats_df,file.path(outdir_all,'LOO_summary_30_models.csv'),row.names=FALSE)

# Verify the originally reported R=10 experiment's difficult replication 4.
check <- stats_df[stats_df$replication==4,]
reference_bic <- c('BF'=-73.86018,'COMBINED-BF'=48.91001,'BATCH-BF'=-74.442847013)
for (alg in algorithms_all) {
  got <- check$BIC[check$algorithm==alg]
  if (length(got)!=1 || abs(got-reference_bic[[alg]])>0.02)
    stop('Replication 4 BIC differs from published experiment (',alg,')')
}

# Representative selection rule, determined before inspecting new plots:
# Rank replications by difference between CombinedFast and BatchFast LOO RMSE.
# Include low, median, high contrasts; 4 is chosen automatically only if ranked accordingly.
wide <- reshape(stats_df[c('replication','algorithm','LOO_RMSE')],
                idvar='replication',timevar='algorithm',direction='wide')
wide$gap <- wide[['LOO_RMSE.COMBINED-BF']] - wide[['LOO_RMSE.BATCH-BF']]
wide <- wide[order(wide$gap,wide$replication),]
representatives <- as.integer(wide$replication[c(1,5,10)])
others <- setdiff(1:10,representatives)
cat('Selected representative replications (low/median/high contrast):',
    paste(representatives,collapse=', '),'\n')
write.csv(data.frame(replication=representatives,
                     contrast=c('low','median','high')),
          file.path(outdir_all,'representative_replications.csv'),row.names=FALSE)

# Uniform y-axis within each replication's three method panels;
# ranges can vary ACROSS replications so normal cases remain legible.
plot_one <- function(rep_id,alg) {
  dd <- points_df[points_df$replication==rep_id & points_df$algorithm==alg,]
  rr <- points_df[points_df$replication==rep_id,]
  ss <- stats_df[stats_df$replication==rep_id & stats_df$algorithm==alg,]
  if (nrow(dd)!=50 || nrow(ss)!=1) stop('Missing results for rep ',rep_id,', ',alg)
  xr <- range(rr$fitted,finite=TRUE); yr <- range(c(0,rr$loo_error),finite=TRUE)
  xr <- extendrange(xr,f=0.05); yr <- extendrange(yr,f=0.09)
  jj <- match(alg,algorithms_all)
  plot(dd$fitted,dd$loo_error,xlim=xr,ylim=yr,pch=pch_all[jj],
       col=colors_all[jj],cex=0.62,
       xlab='Fitted value',ylab='LOO residual',
       main=labels_all[jj],cex.main=0.93,cex.axis=0.76,
       cex.lab=0.83)
  abline(h=0,lty=2,col='gray55',lwd=0.7)
  mtext(sprintf('Rep %d  |  BIC %.2f  |  LOO RMSE %.3f',
                rep_id,ss$BIC,ss$LOO_RMSE),side=3,line=0.15,
        cex=0.57,adj=0)
}

plot_grid <- function(reps,file,width,height) {
  grDevices::pdf(file,width=width,height=height,pointsize=9)
  old <- par(no.readonly=TRUE)
  on.exit({par(old);grDevices::dev.off()},add=TRUE)
  par(mfrow=c(length(reps),3),mar=c(3.4,3.6,2.5,0.9),
      oma=c(0.6,0.7,0.5,0.4),mgp=c(2.1,0.7,0))
  for(rep_id in reps) for(alg in algorithms_all) plot_one(rep_id,alg)
}

for (rep_id in 1:10) {
  plot_grid(rep_id,file.path(outdir_all,sprintf('rep%02d_3panels.pdf',rep_id)),
            width=12.5,height=3.9)
}

plot_grid(representatives,file.path(outdir_all,'Figure_Main_3x3.pdf'),
          width=12.5,height=9.2)
plot_grid(others,file.path(outdir_all,'Figure_Appendix_7x3_full.pdf'),
          width=12.5,height=20.8)
plot_grid(others[1:4],file.path(outdir_all,'Figure_Appendix_part1_4x3.pdf'),
          width=12.5,height=12.4)
plot_grid(others[5:7],file.path(outdir_all,'Figure_Appendix_part2_3x3.pdf'),
          width=12.5,height=9.2)
cat('All done. Figures and CSVs saved in: ', normalizePath(outdir_all), '\n')
