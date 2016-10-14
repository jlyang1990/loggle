# Grid search function for LGGM ###############################################################################################
###############################################################################################################################

# Input ###
# X: a p by N matrix containing list of observations
# pos: position of time points where graphs are estimated
# h: bandwidth in kernel function used to generate correlation matrices
# d.list: list of widths of neighborhood
# lambda.list: list of tuning parameters of Lasso penalty
# fit.type: 0: graphical Lasso estimation, 1: pseudo likelihood estimation, 2: sparse partial correlation estimation
# refit.type: 0: likelihood estimation, 1: pseudo likelihood estimation
# early.stop.thres: grid search stops when number of detected edges exceeds early.stop.thres times number of nodes
# epi.abs: list of absolute tolerances in ADMM stopping criterion
# epi.rel: list of relative tolerances in ADMM stopping criterion
# detrend: whether to detrend each variable in data matrix by subtracting kernel weighted moving average
# fit.corr: whether to use sample correlation matrix rather than sample covariance matrix in model fitting
# num.thread: number of threads

# Output ###
# cv.result.list: results from LGGM.combine.cv

LGGM.grid <- function(X, pos = 1:ncol(X), h = 0.8*ncol(X)^(-1/5), 
                      d.list = c(0, 0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 1), 
                      lambda.list = seq(0.15, 0.35, length = 11), fit.type = "glasso", refit.type = "likelihood", 
                      early.stop.thres = 5, epi.abs = ifelse(nrow(X) >= 400, 1e-4, 1e-5), 
                      epi.rel = ifelse(nrow(X) >= 400, 1e-2, 1e-3), detrend = TRUE, fit.corr = TRUE, num.thread = 1) {
  
  p <- dim(X)[1]
  N <- dim(X)[2]
  K <- length(pos)
  D <- length(d.list)
  L <- length(lambda.list)
  
  if(fit.type == "glasso") {
    fit.type <- 0
  } else if(fit.type == "pseudo") {
    fit.type <- 1
  } else if(fit.type == "space") {
    fit.type <- 2
  } else {
    stop("fit.type must be 'glasso', 'pseudo' or 'space'!")
  }
  
  if(refit.type == "likelihood") {
    refit.type <- 0
  } else if(refit.type == "pseudo") {
    refit.type <- 1
  } else {
    stop("refit.type must be 'likelihood' or 'pseudo'!")
  }
  
  if(any(!pos %in% 1:N)) {
    stop("pos must be a subset of 1, 2, ..., N!")
  }
  if(length(h) != 1) {
    stop("h must be a scalar!")
  }
  
  d.list.global <- d.list >= 0.5
  if(any(d.list.global)) {
    d.list <- c(d.list[!d.list.global], 1)
  }
  d.list <- sort(d.list)
  cat("Using d.list:", d.list, "\n")
  
  lambda.list <- sort(lambda.list)
  cat("Using lambda.list:", lambda.list, "\n")
  
  if(length(epi.abs) == 1) {
    epi.abs <- rep(epi.abs, D)
  }
  if(length(epi.rel) == 1) {
    epi.rel <- rep(epi.rel, D)
  }
  
  if(detrend) {
    cat("Detrending each variable in data matrix...\n")
    X <- dataDetrend(X)
  }
    
  cat("Estimating graphs...\n")
    
  cv.result <- LGGM.combine.cv(X, 1:N, pos, h, d.list, lambda.list, fit.type, refit.type, early.stop.thres, epi.abs, epi.rel,
                               fit.corr, cv.thread)
  
  return(cv.result)
}


# optimal model selection (F1 score) #########################################################################################################################################################################################################################

optimal.select = function(pos, d.l, lambda.c, d.pos, edge, overlap.edge, edge.t){
  
  L = dim(edge)[1]; D = length(d.pos); K = dim(edge)[3]; N = length(edge.t)
  
  d.l = d.l[d.pos]; edge = edge[, d.pos, , drop=F]; overlap.edge = overlap.edge[, d.pos, , drop=F]
  
  FDR = array(0, c(L, D, K)); power = array(0, c(L, D, K))
  for(j in 1:D){for(l in 1:L){
    FDR[l, j, ] = 1-overlap.edge[l, j, ]/edge[l, j, ]; power[l, j, ] = overlap.edge[l, j, ]/edge.t[pos]
  }}
  FDR[is.nan(FDR)] = 1
  
  F1 = 2*(1-FDR)*power/(1-FDR+power); F1[is.nan(F1)] = 0
  
  lambda_d.opt = matrix(0, 2, K); d.opt = rep(0, K+1)
  FDR.opt = matrix(0, 3, 2); power.opt = matrix(0, 3, 2)
  
  lambda_d.index = matrix(0, 2, K)
  for(k in 1:K){
    index = which(matrix(F1[, , k], L, D) == max(F1[, , k]), arr.ind = T)
    if(nrow(index)>1){index = matrix(index[nrow(index), ], c(1, 2))}
    lambda_d.index[, k] = index
    lambda_d.opt[, k] = c(d.l[index[2]], lambda.c[index[1]])
  }
  
  d.index = which.max(sapply(1:D, function(d) sum(apply(F1[, d, ], 2, max))))
  d.lambda.index = sapply(1:K, function(k) which.max(F1[, d.index, k]))
  d.opt[1] = d.l[d.index]; d.opt[-1] = lambda.c[d.lambda.index]
  
  F1.avg = apply(F1, c(1, 2), mean)
  lambda_d.avg.index = which(F1.avg == max(F1.avg), arr.ind = T)
  lambda_d.avg.opt = c(d.l[lambda_d.avg.index[2]], lambda.c[lambda_d.avg.index[1]])
  
  FDR.opt[, 1] = c(mean(sapply(1:K, function(k) FDR[lambda_d.index[1, k], lambda_d.index[2, k], k])), mean(sapply(1:K, function(k) FDR[d.lambda.index[k], d.index, k])), mean(FDR[lambda_d.avg.index[1], lambda_d.avg.index[2], ]))
  FDR.opt[, 2] = c(sd(sapply(1:K, function(k) FDR[lambda_d.index[1, k], lambda_d.index[2, k], k])),sd(sapply(1:K, function(k) FDR[d.lambda.index[k], d.index, k])), sd(FDR[lambda_d.avg.index[1], lambda_d.avg.index[2], ]))/sqrt(K)
  power.opt[, 1] = c(mean(sapply(1:K, function(k) power[lambda_d.index[1, k], lambda_d.index[2, k], k])), mean(sapply(1:K, function(k) power[d.lambda.index[k], d.index, k])), mean(power[lambda_d.avg.index[1], lambda_d.avg.index[2], ]))
  power.opt[, 2] = c(sd(sapply(1:K, function(k) power[lambda_d.index[1, k], lambda_d.index[2, k], k])),sd(sapply(1:K, function(k) power[d.lambda.index[k], d.index, k])), sd(power[lambda_d.avg.index[1], lambda_d.avg.index[2], ]))/sqrt(K)
  F1.opt = c(mean(sapply(1:K, function(k) F1[lambda_d.index[1, k], lambda_d.index[2, k], k])), mean(sapply(1:K, function(k) F1[d.lambda.index[k], d.index, k])), mean(F1[lambda_d.avg.index[1], lambda_d.avg.index[2], ]))
  
  return(list(lambda_d.opt, d.opt, lambda_d.avg.opt, FDR.opt, power.opt, F1.opt))
}


# summary with underlying truth #########################################################################################################################################################################################################################