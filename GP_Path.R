# DUAL-PATH FAST BAS: beta-first vs theta-first, choose smaller full-profile BIC
# Same functions, same LHS samples, same GP profile ML, same jitter.
# BF fully fits each theta neighbor; COMBINED screens beta and theta; BATCH also tries multi-variable changes.
# Reference BF uses beta top 2; COMBINED/BATCH use the configurable shortcuts.
# Based on BAS_new_tests_R5.R. Base R only. Pilot: no global optimality claim.
# source("GP_Path.R")
# Outputs go to a NEW folder; previous study files are untouched.

## SETTINGS: targeted difficult replicates first; no n=100 in this project
R_pilot <- 10L
rep_ids_by_case <- list(baseline4=c(1L), sobol6=c(2L,7L), sparse_sobol10=c(4L))
n_grid <- c(50L)
case_names <- c('baseline4','sobol6','sparse_sobol10')
main_jitter <- 1e-10
ntest <- 300L
outdir <- 'BAS_DualPath_Diagnostic_results'
seed0 <- 20261106L # identical training/test design to prior R=1 pilot
# The original beta screen retained for BOTH methods in this experiment.
fast_topk <- 1L
# Reference beta candidate shortlist is held constant across both algorithms.
# Set TRUE to test every beta neighbor, much slower.
reference_exhaustive <- FALSE
reference_topk <- 2L
# Optimization configuration matches old pilot.
nstarts <- 1L
maxit <- 90L
max_moves <- 40L
# Independent runs per replicate; time includes fresh cache and GP start.
algorithms <- c('BATCH-BF','DUAL-BATCH') # optional: 'DUAL-COMBINED', 'BF'
# Optional exact-GLS algebra self-test before Monte Carlo.
check_gls_formula <- TRUE
# Theta screening: rank neighbors by profiled GLS-BIC with inherited ell;
# full covariance MLE for top K only. 2 is conservative pilot.
theta_topk <- 2L
# Test additional theta lengths for added coordinates.
theta_add_ell <- c(0.2, 1.0, 5.0)
# Reference BF always fits every theta neighbor (as in original algorithm).
theta_reference_exhaustive <- TRUE
# Combined: beta and theta use top K fixed-covariance screenings.
combined_beta_topk <- 2L
# Batch: screen a small candidate pool, permit 2 or 3 variables per move.
batch_pool <- 4L
batch_max_size <- 3L
batch_beta_topk <- 2L
batch_theta_topk <- 2L

stopifnot(R_pilot>=1L,all(n_grid>=3L),reference_topk>=1L,theta_topk>=1L,
          batch_pool>=2L,batch_max_size>=2L,batch_beta_topk>=1L,batch_theta_topk>=1L)
dir.create(outdir,showWarnings=FALSE,recursive=TRUE)

random_lhs_base <- function(n,d) {
  M <- matrix(NA_real_,n,d)
  for(j in seq_len(d)) M[,j] <- (sample.int(n)-runif(n))/n
  as.data.frame(M)
}
# Fixed inputs on [0,1]^d for every example.
evaluate_test <- function(D,case) {
  XX <- as.matrix(D)
  d <- ncol(XX)
  if(case=='baseline4') return(1 + 1.5*XX[,1] -0.7*XX[,2] +
                               0.4*sin(2*pi*XX[,3])+0.25*XX[,1]*XX[,3])
  if(case=='sobol6') {
    a <- seq_len(d) # a_j=j, dimension 6 (Marrel 2008 setting)
    return(apply(sweep(abs(4*XX-2),2,a,'+'),1,function(z) prod(z/(1+a))))
  }
  if(case=='sparse_sobol10') {
    a <- c(1,1,rep(500,d-2L))
    return(apply(sweep(abs(4*XX-2),2,a,'+'),1,function(z) prod(z/(1+a))))
  }
  if(case=='kim_product10') {
    a <- c(rep(1,4L),rep(500,d-4L))
    return(-apply(sweep(4*XX-2,2,a,'+'),1,function(z) prod(z/(1+a))))
  }
  stop('unknown test function: ',case)
}
case_dim <- function(case) switch(case,baseline4=4L,sobol6=6L,
                                  sparse_sobol10=10L,kim_product10=10L,
                                  stop('unknown case'))
make_sample <- function(case,n,rep_id) {
  # Common random seeds across methods and jitters for fair pairing.
  ci <- match(case,case_names)
  set.seed(seed0+10000L*ci+100L*n+rep_id)
  d <- case_dim(case)
  X <- random_lhs_base(n,d); Xtest <- random_lhs_base(ntest,d)
  names(X) <- names(Xtest) <- paste0('x',seq_len(d))
  list(X=X,Xtest=Xtest,y=evaluate_test(X,case),ytest=evaluate_test(Xtest,case))
}

cfg <- list(n=50L, ntest=ntest,d=4L,jitter=main_jitter,
            ell_lower=0.003,ell_upper=20,nstarts=nstarts,maxit=maxit,
            max_moves=max_moves,exhaustive_neighbors=FALSE,
            shortlist=reference_topk,tolerance=1e-6,
            selection_mode='reference',reference_exhaustive=reference_exhaustive,
            reference_topk=reference_topk,fast_topk=fast_topk,
            theta_topk=theta_topk,theta_reference_exhaustive=theta_reference_exhaustive,
            combined_beta_topk=combined_beta_topk,batch_pool=batch_pool,
            batch_max_size=batch_max_size,batch_beta_topk=batch_beta_topk,
            batch_theta_topk=batch_theta_topk)
keyset <- function(z) paste(sort(z),collapse=',')
keymodel <- function(b,t) paste0('B[',keyset(b),']T[',keyset(t),']')
make_F <- function(D,b) {
  cbind('(Intercept)'=rep(1,nrow(D)),
        if(length(b)) as.matrix(D[,b,drop=FALSE]))
}
cov_matrix <- function(A,B,t,ell) {
  AA <- as.matrix(A[,t,drop=FALSE]); BB <- as.matrix(B[,t,drop=FALSE])
  H <- matrix(0,nrow(AA),nrow(BB))
  for(j in seq_along(t)) H <- H +
    outer(AA[,j],BB[,j],'-')^2/(2*ell[j]^2)
  exp(-H)
}

# Profile full ML (not REML). Fixed, dimensionless diagonal jitter added to R.
profile_gp <- function(ell,b,t) {
  if(length(t)==0 || length(ell)!=length(t) || any(!is.finite(ell)) ||
     any(ell<=0)) return(NULL)
  F <- make_F(X,b)
  R <- cov_matrix(X,X,t,ell)
  diag(R) <- diag(R) + cfg$jitter
  U <- tryCatch(chol(R),error=function(e) NULL)
  if(is.null(U)) return(NULL)
  Ri <- function(M) backsolve(U,forwardsolve(t(U),as.matrix(M)))
  RF <- Ri(F); Ry <- Ri(y)
  A <- crossprod(F,RF)
  beta <- tryCatch(as.vector(solve(A,crossprod(F,Ry))),error=function(e) NULL)
  if(is.null(beta) || any(!is.finite(beta))) return(NULL)
  res <- as.vector(y-F%*%beta)
  alpha <- as.vector(Ri(res))
  sig2 <- sum(res*alpha)/length(y)
  if(!is.finite(sig2) || sig2<=0) return(NULL)
  ll <- -.5*(length(y)*(log(2*pi)+1+log(sig2))+
              2*sum(log(diag(U))))
  if(!is.finite(ll)) return(NULL)
  list(b=b,t=t,ell=ell,theta=setNames(1/(2*ell^2),t),
       beta=beta,sigma2=sig2,logLik=ll,alpha=alpha,
       R=R,F=F,chol=U)
}

predict_gp <- function(fit,Dnew) {
  crosscov <- cov_matrix(Dnew,X,fit$t,fit$ell)
  as.vector(make_F(Dnew,fit$b)%*%fit$beta + crosscov%*%fit$alpha)
}

# Models are cached to avoid fitting the same subset repeatedly.
cache <- new.env(parent=emptyenv())
fit_calls <- 0L
fit_model <- function(b,t,parent=NULL,refresh=FALSE) {
  b <- sort(unique(b)); t <- sort(unique(t))
  if(!length(t)) return(NULL) # Gaussian GP with zero nugget cannot have t={}
  key <- keymodel(b,t)
  if(!refresh && exists(key,envir=cache,inherits=FALSE))
    return(get(key,envir=cache,inherits=FALSE))
  p <- length(t)
  lo <- cfg$ell_lower; hi <- cfg$ell_upper
  objective <- function(logell) {
    z <- tryCatch(profile_gp(exp(logell),b,t),error=function(e) NULL)
    if(is.null(z)) return(1e15)
    -z$logLik
  }
  grid <- c(0.01,0.04,0.15,0.5,1.5,5)
  starts <- lapply(grid,function(v) rep(v,p))
  if(!is.null(parent)) {
    parent_ell <- setNames(parent$ell,parent$t)
    st <- vapply(t,function(j) {
      if(j %in% names(parent_ell)) parent_ell[[j]] else 0.25
    },numeric(1))
    starts <- c(list(st),starts)
  }
  if(p>1) {
    st <- rep(0.2,p); st[1] <- 1.5
    starts <- c(starts,list(st))
  }
  starts <- lapply(starts,function(z) pmax(lo,pmin(hi,z)))
  scores <- vapply(starts,function(z) objective(log(z)),numeric(1))
  viable <- which(is.finite(scores) & scores < 1e14)
  if(!length(viable)) {
    assign(key,NULL,envir=cache)
    return(NULL)
  }
  viable <- viable[order(scores[viable])]
  # Always test warm start plus the best available grid starts.
  use <- viable[seq_len(min(length(viable),cfg$nstarts+as.integer(!is.null(parent))))]
  best <- NULL
  for(ix in use) {
    op <- tryCatch(optim(log(starts[[ix]]),objective,method='L-BFGS-B',
       lower=rep(log(lo),p),upper=rep(log(hi),p),
       control=list(maxit=cfg$maxit,factr=1e7)),error=function(e) NULL)
    if(is.null(op) || !is.finite(op$value) || op$value>=1e14) next
    cand <- tryCatch(profile_gp(exp(op$par),b,t),error=function(e) NULL)
    if(is.null(cand)) next
    cand$convergence <- op$convergence
    if(is.null(best) || cand$logLik>best$logLik) best <- cand
  }
  if(is.null(best)) {
    best <- profile_gp(starts[[viable[1]]],b,t)
    if(!is.null(best)) best$convergence <- 99L
  }
  if(!is.null(best)) {
    best$k <- 1L+length(b)+length(t)+1L # intercept + beta + theta + sigma2
    best$BIC <- -2*best$logLik+best$k*log(cfg$n)
    best$rcond <- tryCatch(rcond(best$R),error=function(e) NA_real_)
    best$boundary <- any(best$ell < lo*1.001 | best$ell > hi*0.999)
    best$key <- key
    fit_calls <<- fit_calls+1L
  }
  assign(key,best,envir=cache)
  best
}

neighbors <- function(fit,block=c('beta','theta','both')) {
  block <- match.arg(block)
  out <- list()
  if(block %in% c('beta','both')) for(j in allvars) {
    b <- if(j %in% fit$b) setdiff(fit$b,j) else union(fit$b,j)
    out[[length(out)+1L]] <- list(b=b,t=fit$t,block='beta',variable=j)
  }
  if(block %in% c('theta','both')) for(j in allvars) {
    t <- if(j %in% fit$t) setdiff(fit$t,j) else union(fit$t,j)
    if(length(t)) out[[length(out)+1L]] <-
      list(b=fit$b,t=t,block='theta',variable=j)
  }
  out
}

# Exact CONDITIONAL Gaussian ML BIC for any beta subset with R and ell fixed.
# One whitening operation per step; reused across all beta neighbors.
# Whitened design has intercept plus all d predictors.
prepare_gls_scores <- function(fit) {
  U <- fit$chol
  Ffull <- cbind(1,as.matrix(X[,allvars,drop=FALSE]))
  colnames(Ffull) <- c('(Intercept)',allvars)
  whitened <- forwardsolve(t(U),cbind(y,Ffull))
  wy <- whitened[,1L]
  wf <- whitened[,-1L,drop=FALSE]
  # forwardsolve() can discard dimnames. Restore the original design labels.
  colnames(wf) <- colnames(Ffull)
  Qbase <- cfg$n*fit$sigma2
  function(bnew) {
    cols <- match(c('(Intercept)',bnew),colnames(wf))
    if(anyNA(cols)) stop('GLS screening: unknown beta variable: ',
                         paste(c('(Intercept)',bnew)[is.na(cols)],collapse=','))
    Z <- wf[,cols,drop=FALSE]
    bh <- tryCatch(qr.coef(qr(Z,tol=1e-12),wy),error=function(e) NULL)
    if(is.null(bh) || any(!is.finite(bh))) return(Inf)
    res <- wy - as.vector(Z%*%bh)
    Qnew <- sum(res^2)
    if(!is.finite(Qnew) || Qnew<=0 || Qbase<=0) return(Inf)
    cfg$n*log(Qnew/Qbase) + (length(bnew)-length(fit$b))*log(cfg$n)
  }
}

# Cheap, model-specific theta proxy: profile beta and sigma2 analytically at
# inherited covariance lengths (no theta optimization). For theta deletion,
# retained lengths come directly from the parent fit. For addition, test
# a small fixed set for the added coordinate; no GP MLE is performed.
# This is a HEURISTIC ranker, NOT the reoptimized BIC.
proxy_theta_bic <- function(cur,m) {
  old <- setNames(cur$ell,cur$t)
  added <- setdiff(m$t,cur$t)
  if(length(added)>1L) stop('theta proxy expects single-variable moves')
  candidates <- if(length(added)) theta_add_ell else NA_real_
  vals <- vapply(candidates,function(v) {
    ell <- vapply(m$t,function(j) {
      if(j %in% names(old)) old[[j]] else v
    },numeric(1))
    z <- tryCatch(profile_gp(ell,m$b,m$t),error=function(e) NULL)
    if(is.null(z)) return(Inf)
    k <- 2L+length(m$b)+length(m$t)
    -2*z$logLik+k*log(cfg$n)
  },numeric(1))
  min(vals)
}

# Candidate selection. All proxy BIC values are conditional on inherited
# covariance lengths, with beta and sigma^2 profiled. A proxy may rank an
# ultimately better model poorly: screening is HEURISTIC.
# In batch mode, all single changes remain eligible alongside limited
# multi-variable changes, including additions/deletions of beta variables
# and simultaneous deletions (not additions) of theta variables.
unique_models <- function(ms) {
  keys <- vapply(ms,function(m) keymodel(m$b,m$t),character(1))
  ms[!duplicated(keys)]
}
make_batch_neighbors <- function(current,block) {
  singles <- neighbors(current,block)
  if(!length(singles)) return(singles)
  if(block=='beta') {
    score <- prepare_gls_scores(current)
    order1 <- order(vapply(singles,function(m) score(m$b),numeric(1)))
    topvars <- vapply(singles[head(order1,cfg$batch_pool)],function(m) m$variable,character(1))
    out <- singles
    if(length(topvars)>=2L && cfg$batch_max_size>=2L) for(k in seq.int(2L,min(cfg$batch_max_size,length(topvars)))) {
      for(v in combn(topvars,k,simplify=FALSE)) {
        bnew <- sort(unique(c(setdiff(current$b,v),setdiff(v,current$b))))
        out[[length(out)+1L]] <- list(b=bnew,t=current$t,block='beta',
                                     variable=paste(v,collapse='+'))
      }
    }
    return(unique_models(out))
  }
  # Theta batch candidates: screen only currently active coordinates.
  removals <- Filter(function(m) length(m$t)<length(current$t),singles)
  if(length(removals)<2L) return(singles)
  ranking <- order(vapply(removals,function(m) proxy_theta_bic(current,m),numeric(1)))
  pool <- vapply(removals[head(ranking,cfg$batch_pool)],function(m) m$variable,character(1))
  out <- singles
  if(min(length(pool),length(current$t)-1L,cfg$batch_max_size)>=2L) for(k in seq.int(2L,min(cfg$batch_max_size,length(pool),length(current$t)-1L))) {
    for(v in combn(pool,k,simplify=FALSE)) {
      out[[length(out)+1L]] <- list(b=current$b,t=setdiff(current$t,v),
                                   block='theta',variable=paste(v,collapse='+'))
    }
  }
  unique_models(out)
}
choose_best_move <- function(current,block) {
  mode <- cfg$selection_mode
  nb <- if(mode=='batch') make_batch_neighbors(current,block) else neighbors(current,block)
  if(!length(nb)) return(NULL)
  ids <- seq_along(nb)
  if(block=='beta') {
    if(mode=='reference' && cfg$reference_exhaustive) {
      ids <- seq_along(nb)
    } else {
      score <- prepare_gls_scores(current)
      delta <- vapply(nb,function(m) score(m$b),numeric(1))
      k <- switch(mode,reference=cfg$reference_topk,
                  combined=cfg$combined_beta_topk,
                  batch=cfg$batch_beta_topk,
                  cfg$reference_topk)
      ids <- head(order(delta),k)
    }
  } else if(block=='theta' && mode %in% c('combined','batch')) {
    proxy <- vapply(nb,function(m) proxy_theta_bic(current,m),numeric(1))
    k <- if(mode=='batch') cfg$batch_theta_topk else cfg$theta_topk
    ids <- head(order(proxy),k)
  }
  # Always compare fully refitted BIC to current (same fixed jitter).
  best <- NULL
  for(i in ids) {
    m <- nb[[i]]
    f <- fit_model(m$b,m$t,parent=current)
    if(is.null(f) || !is.finite(f$BIC)) next
    if(f$BIC < current$BIC-cfg$tolerance &&
       (is.null(best) || f$BIC<best$BIC)) best <- f
  }
  best
}

trace_rows <- list()
append_trace <- function(method,move,fit) {
  trace_rows[[length(trace_rows)+1L]] <<- data.frame(
    method=method,step=move,key=fit$key,BIC=fit$BIC,
    Sbeta=keyset(fit$b),Stheta=keyset(fit$t),
    logLik=fit$logLik,convergence=fit$convergence,
    rcond=fit$rcond,boundary=fit$boundary)
}
block_search <- function(cur,block,method,steps_left) {
  used <- 0L
  repeat {
    if(used>=steps_left) break
    next_fit <- choose_best_move(cur,block)
    if(is.null(next_fit)) break
    cur <- next_fit; used <- used+1L
    append_trace(method,paste0(block,'-',used),cur)
  }
  list(fit=cur,moves=used)
}
run_search <- function(method,start) {
  # BF or TF denotes a SINGLE pass through the two selection blocks;
  # both start from the SAME initial model (intercept, all theta variables).
  now <- start; step <- 0L; append_trace(method,'start',now)
  order_blocks <- if(method=='BF') c('beta','theta') else
                  if(method=='TF') c('theta','beta') else stop('Unknown search order')
  for(bl in order_blocks) {
    if(step>=cfg$max_moves) break
    z <- block_search(now,bl,method,cfg$max_moves-step)
    now <- z$fit; step <- step+z$moves
  }
  list(fit=now,moves=step)
}

# Small deterministic validation of GLS conditional objective on a common R.
verify_gls <- function() {
  dat <- make_sample('sobol6',50L,1L)
  cfg$n <<- 50L; cfg$d <<- 6L; cfg$jitter <<- main_jitter
  X <<- dat$X; y <<- dat$y
  allvars <<- paste0('x',1:6)
  b <- c('x1','x3'); t <- c('x1','x2','x4')
  ell <- c(0.25,0.55,0.3)
  f <- profile_gp(ell,b,t)
  if(is.null(f)) stop('GLS self-test: initial profile invalid')
  a <- prepare_gls_scores(f)
  cand <- list(c('x1','x3','x5'),c('x1'),c('x2','x3'),c('x1','x2','x3','x5'))
  for(bnew in cand) {
    direct <- profile_gp(ell,bnew,t)
    if(is.null(direct)) stop('GLS self-test: direct profile invalid')
    kdiff <- length(bnew)-length(b)
    true_delta <- -2*(direct$logLik-f$logLik)+kdiff*log(50)
    approx <- a(bnew)
    if(!is.finite(approx) || abs(approx-true_delta)>1e-5)
      stop('GLS test failed: expected=',true_delta,' got=',approx)
  }
  cat('GLS conditional BIC self-check: PASS (4 beta moves, including batch)\n')
}
if(check_gls_formula) verify_gls()

# ------------------- DUAL PATH DIAGNOSTIC DRIVER -----------------------
# Distinct outputs; previously completed study directories are untouched.
# Example to also diagnose sparse_sobol10 rep 2/3/6/10:
# rep_ids_by_case$sparse_sobol10 <- c(2L,3L,4L,6L,10L)
# Each of the two paths has a fresh cache; DUAL runtime = SUM of both paths.
# For a same-replicate BF reference WITHOUT refitting BF, supply its individual.csv.
reference_csv <- ''
resume_folder <- ''

stopifnot(all(algorithms %in% c('BATCH-BF','DUAL-BATCH','DUAL-COMBINED','BF')))
if(nzchar(resume_folder)) {
  if(!dir.exists(resume_folder)) stop('resume_folder does not exist: ',resume_folder)
  runfolder <- normalizePath(resume_folder,winslash='/',mustWork=TRUE)
} else {
  runfolder <- file.path(outdir,paste0('run_',format(Sys.time(),'%Y%m%d_%H%M%S')))
  if(dir.exists(runfolder)) stop('Output folder exists: ',runfolder)
  dir.create(runfolder,recursive=TRUE)
}
individual_file <- file.path(runfolder,'individual.csv')
paths_file <- file.path(runfolder,'path_results.csv')
trace_file <- file.path(runfolder,'traces.csv')
write_one <- function(z,filename) {
  if(is.null(z) || !nrow(z)) return(invisible(NULL))
  write.table(z,filename,sep=',',row.names=FALSE,quote=TRUE,
              col.names=!file.exists(filename),append=file.exists(filename),na='')
}
key_run <- function(case,n,r,alg) paste(case,n,r,alg,sep='|')
read_completed <- function() {
  if(!file.exists(individual_file)) return(character(0))
  z <- tryCatch(read.csv(individual_file,stringsAsFactors=FALSE),error=function(e) NULL)
  if(is.null(z) || !nrow(z)) return(character(0))
  unique(mapply(key_run,z$case,z$n,z$replicate,z$algorithm,USE.NAMES=FALSE))
}

# One path, fresh cache, same GP fitting as original R10 study.
fit_path <- function(dat,case,n,rep_id,mode,path,algorithm) {
  cfg$n <<- n; cfg$d <<- case_dim(case); cfg$jitter <<- main_jitter
  cfg$selection_mode <<- mode
  allvars <<- paste0('x',seq_len(cfg$d))
  X <<- dat$X; y <<- dat$y; Xtest <<- dat$Xtest; ytest <<- dat$ytest
  cache <<- new.env(parent=emptyenv()); fit_calls <<- 0L
  trace_rows <<- list()
  t0 <- proc.time()[['elapsed']]
  err <- ''
  result <- tryCatch({
    initial <- fit_model(character(0),allvars)
    if(is.null(initial)) stop('initial GP fit failed')
    z <- run_search(path,initial)
    f <- z$fit
    list(f=f,moves=z$moves,rmspe=sqrt(mean((predict_gp(f,Xtest)-ytest)^2)))
  },error=function(e) {err <<- conditionMessage(e); NULL})
  seconds <- proc.time()[['elapsed']]-t0
  if(is.null(result)) {
    row <- data.frame(case=case,n=n,replicate=rep_id,algorithm=algorithm,path=path,
        Sbeta=NA_character_,Stheta=NA_character_,BIC=NA_real_,RMSPE=NA_real_,
        seconds=seconds,number_fits=fit_calls,moves=NA_integer_,
        convergence=NA_integer_,rcond=NA_real_,boundary=NA,warning=err)
    tr <- NULL
  } else {
    f <- result$f
    row <- data.frame(case=case,n=n,replicate=rep_id,algorithm=algorithm,path=path,
        Sbeta=keyset(f$b),Stheta=keyset(f$t),BIC=f$BIC,RMSPE=result$rmspe,
        seconds=seconds,number_fits=fit_calls,moves=result$moves,
        convergence=f$convergence,rcond=f$rcond,boundary=f$boundary,warning='')
    tr <- do.call(rbind,trace_rows)
    tr$case <- case; tr$n <- n; tr$replicate <- rep_id
    tr$algorithm <- algorithm; tr$path <- path; tr$jitter <- main_jitter
  }
  list(row=row,trace=tr)
}

run_one <- function(case,n,rep_id,algorithm,dat) {
  mode <- if(algorithm %in% c('BATCH-BF','DUAL-BATCH')) 'batch' else
          if(algorithm=='DUAL-COMBINED') 'combined' else 'reference'
  paths <- if(startsWith(algorithm,'DUAL-')) c('BF','TF') else 'BF'
  results <- lapply(paths,function(path)
    fit_path(dat,case,n,rep_id,mode,path,algorithm))
  pr <- do.call(rbind,lapply(results,function(z) z$row))
  write_one(pr,paths_file)
  for(z in results) write_one(z$trace,trace_file)
  good <- which(is.finite(pr$BIC) & is.finite(pr$RMSPE))
  if(length(good)) {
    best <- good[which.min(pr$BIC[good])]
    chosen <- pr[best,,drop=FALSE]
  } else {
    chosen <- pr[1L,,drop=FALSE]
    chosen$path <- NA_character_
    chosen$BIC <- NA_real_; chosen$RMSPE <- NA_real_
  }
  # Runtime and number_fits include BOTH independent paths for a dual method.
  chosen$seconds <- sum(pr$seconds)
  chosen$number_fits <- sum(pr$number_fits)
  chosen$path_BF_BIC <- if('BF' %in% paths) pr$BIC[which(paths=='BF')] else NA_real_
  chosen$path_TF_BIC <- if('TF' %in% paths) pr$BIC[which(paths=='TF')] else NA_real_
  chosen$BIC_gain_TF_over_BF <- if(length(paths)==2L && all(is.finite(pr$BIC)))
        pr$BIC[1L]-min(pr$BIC) else NA_real_
  chosen$jitter <- main_jitter
  write_one(chosen,individual_file)
  cat(sprintf('%-16s rep=%2d %-14s winner=%-2s BIC=%9.2f RMSPE=%9.5g time=%6.1fs fits=%d | BFpath=%9.2f TFpath=%9.2f\n',
    case,rep_id,algorithm,ifelse(is.na(chosen$path),'--',chosen$path),chosen$BIC,
    chosen$RMSPE,chosen$seconds,chosen$number_fits,
    chosen$path_BF_BIC,chosen$path_TF_BIC))
  flush.console()
  invisible(chosen)
}

summarize_now <- function() {
  if(!file.exists(individual_file)) return(invisible(NULL))
  z <- read.csv(individual_file,stringsAsFactors=FALSE)
  if(!nrow(z)) return(invisible(NULL))
  keys <- mapply(key_run,z$case,z$n,z$replicate,z$algorithm,USE.NAMES=FALSE)
  z <- z[!duplicated(keys,fromLast=TRUE),,drop=FALSE]
  write.csv(z,file.path(runfolder,'all_results.csv'),row.names=FALSE,na='')
  sm <- do.call(rbind,lapply(split(z,list(z$case,z$algorithm),drop=TRUE),function(a) {
    ok <- is.finite(a$BIC)&is.finite(a$RMSPE)
    data.frame(case=a$case[1L],algorithm=a$algorithm[1L],R=nrow(a),valid=sum(ok),
      mean_BIC=if(any(ok)) mean(a$BIC[ok]) else NA_real_,
      mean_RMSPE=if(any(ok)) mean(a$RMSPE[ok]) else NA_real_,
      mean_seconds=mean(a$seconds),mean_fits=mean(a$number_fits),
      TF_wins=sum(a$path=='TF',na.rm=TRUE),
      nonzero_convergence=sum(!is.na(a$convergence)&a$convergence!=0L),
      min_rcond=if(any(is.finite(a$rcond))) min(a$rcond,na.rm=TRUE) else NA_real_)
  }))
  rownames(sm) <- NULL
  write.csv(sm,file.path(runfolder,'summary.csv'),row.names=FALSE,na='')
  comps <- list()
  for(case in case_names) for(n in n_grid) {
    a <- z[z$case==case & z$n==n,,drop=FALSE]
    base <- a[a$algorithm=='BATCH-BF',,drop=FALSE]
    dual <- a[a$algorithm=='DUAL-BATCH',,drop=FALSE]
    m <- merge(base,dual,by='replicate',suffixes=c('_single','_dual'))
    if(!nrow(m)) next
    good <- is.finite(m$BIC_single)&is.finite(m$BIC_dual)
    m <- m[good,,drop=FALSE]
    if(!nrow(m)) next
    comps[[length(comps)+1L]] <- data.frame(case=case,n=n,R=nrow(m),
      dual_BIC_better=sum(m$BIC_dual<m$BIC_single-1e-5),
      dual_BIC_same=sum(abs(m$BIC_dual-m$BIC_single)<1e-5),
      mean_BIC_dual_minus_single=mean(m$BIC_dual-m$BIC_single),
      mean_RMSPE_dual_minus_single=mean(m$RMSPE_dual-m$RMSPE_single),
      median_time_ratio_dual_to_single=median(m$seconds_dual/m$seconds_single),
      same_model=mean(m$Sbeta_single==m$Sbeta_dual & m$Stheta_single==m$Stheta_dual))
  }
  if(length(comps)) write.csv(do.call(rbind,comps),file.path(runfolder,'paired.csv'),row.names=FALSE,na='')
  if(nzchar(reference_csv)) {
    if(!file.exists(reference_csv)) warning('reference_csv does not exist: ',reference_csv)
    else {
      ref <- read.csv(reference_csv,stringsAsFactors=FALSE)
      if(all(c('case','n','replicate','algorithm','BIC','RMSPE') %in% names(ref))) {
        ref <- ref[ref$algorithm=='BF',,drop=FALSE]
        m <- merge(z,ref,by=c('case','n','replicate'),suffixes=c('_dual','_oldBF'))
        if(nrow(m)) {
          m$delta_BIC <- m$BIC_dual-m$BIC_oldBF
          m$delta_RMSPE <- m$RMSPE_dual-m$RMSPE_oldBF
          write.csv(m,file.path(runfolder,'vs_old_BF.csv'),row.names=FALSE,na='')
        }
      } else warning('reference_csv missing required columns')
    }
  }
  invisible(NULL)
}
cat('GLS check was completed; comparing beta-first vs theta-first routes.\n')
cat('=== Dual Path Fast BAS diagnostic ===\n')
cat('OUTPUT FOLDER: ',normalizePath(runfolder,winslash='/',mustWork=TRUE),'\n')
cat('TO RESUME: set resume_folder to this exact path and source again.\n')
completed <- read_completed()
for(case in case_names) for(n in n_grid) for(rep_id in rep_ids_by_case[[case]]) {
  dat <- make_sample(case,n,rep_id)
  for(algorithm in algorithms) {
    key <- key_run(case,n,rep_id,algorithm)
    if(key %in% completed) next
    run_one(case,n,rep_id,algorithm,dat)
    completed <- c(completed,key)
    summarize_now()
  }
}
summarize_now()
cat('Completed or already checkpointed. Results: ',normalizePath(runfolder,winslash='/',mustWork=TRUE),'\n')
cat('NOTE: Dual BIC cannot exceed BF-path BIC within the same screening/optimizer implementation,\n')
cat('      but dual may be worse than OLD BF with a different screening configuration.\n')
cat('NOTE: Both paths use their own fresh cache; dual timing includes both paths.\n')
cat('WARNING: top-k and batch pool remain heuristic, no guarantee of best model.\n')
