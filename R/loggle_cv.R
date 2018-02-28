
# Cross validation function for loggle (with h fixed) ##################################################################
########################################################################################################################

# Input ###
# X: a p by N data matrix containing observations on a time grid ranging from 0 to 1
# pos: indices of time points where graphs are estimated
# h: bandwidth in kernel smoothed sample covariance/correlation matrix 
# d.list: a grid of widths of neighborhood centered at each time point specified by pos
# lambda.list: a grid of tuning parameters of lasso penalty at each time point specified by pos
# cv.fold: number of cross-validation folds
# fit.type: likelihood: likelihood estimation, 
#           pseudo: pseudo likelihood estimation, 
#           space: sparse partial correlation estimation
# return.select: if TRUE, return model selection result
# select.type: "all_flexible":optimal d and lambda can vary across time points specified by pos, 
#              "d_fixed": optimal d is fixed and optimal lambda can vary across time points specified by pos, 
#              "all_fixed": optimal d and lambda are fixed across time points specified by pos
# cv.vote.thres: an edge is kept after cv.vote if and only if it exists in no less than cv.vote.thres*cv.fold cv folds
# early.stop.thres: grid search stops when the ratio between edge number and variable number exceeds early.stop.thres
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion
# max.step: maximum steps in ADMM iteration
# detrend: if TRUE, subtract kernel weighted moving average for each variable in data matrix,
#          if FALSE, subtract overall average for each variable in data matrix
# fit.corr: if TRUE, use sample correlation matrix in model fitting,
#           if FALSE, use sample covariance matrix in model fitting
# h.correct: if TRUE, apply bandwidth adjustment for validation sets due to sample size difference
# num.thread: number of threads used in parallel computing
# print.detail: if TRUE, print details in model fitting procedure

# Output ###
# cv.score: an array of cv scores for each combination of d, lambda, time point and cv fold
# cv.result.fold: a list of model fitting results for each cv fold, including:
#                 a list of estimated precision matrices for each combination of d, lambda and time point, 
#                 an array of numbers of edges for each combination of d, lambda and time point, 
#                 a list of edges for each combination of d, lambda and time point
# cv.select.result: results from loggle.cv.select.h if return.select = TRUE

loggle.cv.h <- function(X, pos = 1:ncol(X), h = 0.8*ncol(X)^(-1/5), 
                        d.list = c(0, 0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 1), 
                        lambda.list = seq(0.15, 0.35, 0.02), cv.fold = 5, fit.type = "pseudo",
                        return.select = TRUE, select.type = "all_flexible", cv.vote.thres = 0.8, early.stop.thres = 5, 
                        epi.abs = 1e-4, epi.rel = 1e-2, max.step = 500, detrend = TRUE, fit.corr = TRUE, 
                        h.correct = TRUE, num.thread = 1, print.detail = TRUE) {
  
  p <- dim(X)[1]
  N <- dim(X)[2]
  K <- length(pos)
  L <- length(lambda.list)
  
  if(fit.type == "likelihood") {
    fit.type <- 0
  } else if(fit.type == "pseudo") {
    fit.type <- 1
  } else if(fit.type == "space") {
    fit.type <- 2
  } else {
    stop("fit.type must be 'likelihood', 'pseudo' or 'space'!")
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
  D <- length(d.list)
  cat("Using d.list:", d.list, "\n")
  
  lambda.list <- sort(lambda.list)
  cat("Using lambda.list:", lambda.list, "\n")
  
  if(return.select && !select.type %in% c("all_flexible", "d_fixed", "all_fixed")) {
    stop("select.type must be 'all_flexible', 'd_fixed' or 'all_fixed'!")
  }
  
  if(length(epi.abs) == 1) {
    epi.abs <- rep(epi.abs, D)
  }
  if(length(epi.rel) == 1) {
    epi.rel <- rep(epi.rel, D)
  }
  
  cat("Detrending each variable in data matrix...\n")
  X <- dataDetrend(X, detrend)
  
  if(h.correct) {
    h.test <- h * (cv.fold-1)^(1/5)
  } else {
    h.test <- h
  }
  
  cv.score <- array(NA, c(L, D, K, cv.fold))
  rownames(cv.score) <- lambda.list
  colnames(cv.score) <- d.list
  cv.result.fold <- vector("list", cv.fold)
  
  if(d.list[1] != 1 && num.thread > 1) {
    registerDoParallel(num.thread)
  }
  
  for(i in 1:cv.fold) {
    
    cat("\nRunning fold", i, "out of", cv.fold, "folds...\n")
    
    pos.test <- seq(i, N, cv.fold)
    pos.train <- (1:N)[-pos.test]
    
    result.i <- loggle.combine.cv(X, pos.train, pos, h, d.list, lambda.list, fit.type, early.stop.thres, epi.abs, 
                                  epi.rel, max.step, fit.corr, num.thread, print.detail)
    cv.result.fold[[i]] <- result.i
    
    cat("Calculating cross-validation scores for testing dataset...\n")
    
    Sigma.test <- makeCorr(X, pos.test, h.test, fit.corr = FALSE)$Corr
    
    for(d in 1:D) {
      for(l in 1:L) {
        for(k in 1:K) {
          Omega <- as.matrix(result.i$Omega[[l, d, k]])
          cv.score[l, d, k, i] <- sum(c(t(Sigma.test[, , pos[k]]))*c(Omega)) - log(det(Omega))
        }
      }
    }
    
    rm(Sigma.test)
  }
  
  cv.result <- list(cv.score = cv.score, cv.result.fold = cv.result.fold)
  
  if(return.select) {
    
    cat(sprintf("\nSelecting models based on %d-fold cross-validation results...\n", cv.fold))
    cv.select.result <- loggle.cv.select.h(cv.result, select.type, cv.vote.thres)
    cv.result$cv.select.result <- cv.select.result
  }
  
  return(cv.result)
}


# Cross validation function for loggle #################################################################################
########################################################################################################################

# Input ###
# X: a p by N data matrix containing observations on a time grid ranging from 0 to 1
# pos: indices of time points where graphs are estimated
# h.list: a grid of bandwidths in kernel smoothed sample covariance/correlation matrix
# d.list: a grid of widths of neighborhood centered at each time point specified by pos
# lambda.list: a grid of tuning parameters of lasso penalty at each time point specified by pos
# cv.fold: number of cross-validation folds
# fit.type: likelihood: likelihood estimation, 
#           pseudo: pseudo likelihood estimation, 
#           space: sparse partial correlation estimation
# return.select: if TRUE, return model selection result
# select.type: "all_flexible":optimal d and lambda can vary across time points specified by pos, 
#              "d_fixed": optimal d is fixed and optimal lambda can vary across time points specified by pos, 
#              "all_fixed": optimal d and lambda are fixed across time points specified by pos
# cv.vote.thres: an edge is kept after cv.vote if and only if it exists in no less than cv.vote.thres*cv.fold cv folds
# early.stop.thres: grid search stops when the ratio between edge number and variable number exceeds early.stop.thres
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion
# max.step: maximum steps in ADMM iteration
# detrend: if TRUE, subtract kernel weighted moving average for each variable in data matrix,
#          if FALSE, subtract overall average for each variable in data matrix
# fit.corr: if TRUE, use sample correlation matrix in model fitting,
#           if FALSE, use sample covariance matrix in model fitting
# h.correct: if TRUE, apply bandwidth adjustment for validation sets due to sample size difference
# num.thread: number of threads used in parallel computing
# print.detail: if TRUE, print details in model fitting procedure

# Output ###
# cv.result.h: a list of model fitting results from loggle.cv.h for each h
# cv.select.result: results from loggle.cv.select if return.select = TRUE

loggle.cv <- function(X, pos = 1:ncol(X), h.list = seq(0.1, 0.3, 0.05), 
                      d.list = c(0, 0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 1), 
                      lambda.list = seq(0.15, 0.35, 0.02), cv.fold = 5, fit.type = "pseudo", 
                      return.select = TRUE, select.type = "all_flexible", cv.vote.thres = 0.8, early.stop.thres = 5,
                      epi.abs = 1e-4, epi.rel = 1e-2, max.step = 500, detrend = TRUE, fit.corr = TRUE, 
                      h.correct = TRUE, num.thread = 1, print.detail = TRUE) {
  
  H <- length(h.list)
  cv.result.h <- vector("list", H)
  
  for(h in 1:H) {
    
    cat("\nRunning h =", h.list[h], "...\n")
    cv.result.h[[h]] <- loggle.cv.h(X, pos, h.list[h], d.list, lambda.list, cv.fold, fit.type, return.select = FALSE, 
                                    select.type, cv.vote.thres, early.stop.thres, epi.abs, epi.rel, max.step, detrend, 
                                    fit.corr, h.correct, num.thread, print.detail)
    cv.result.h[[h]]$h <- h.list[h]
  }
  
  cv.result <- list(cv.result.h = cv.result.h)
  
  if(return.select) {
    
    cat(sprintf("\nSelecting models based on %d-fold cross-validation results...\n", cv.fold))
    cv.select.result <- loggle.cv.select(cv.result, select.type, cv.vote.thres)
    cv.result$cv.select.result <- cv.select.result
  }
  
  return(cv.result)
}


# Selection function for cross validation result (with h fixed) ########################################################
########################################################################################################################

# Input ###
# cv.result: results from loggle.cv.h
# select.type: "all_flexible":optimal d and lambda can vary across time points specified by pos, 
#              "d_fixed": optimal d is fixed and optimal lambda can vary across time points specified by pos, 
#              "all_fixed": optimal d and lambda are fixed across time points specified by pos
# cv.vote.thres: an edge is kept after cv.vote if and only if it exists in no less than cv.vote.thres*cv.fold cv folds

# Output ###
# d.opt: a vector of optimal values of d for each estimated graph
# lambda.opt: a vector of optimal values of lambda for each estimated graph
# cv.score.opt: optimal cv score (averaged over time points and cv folds)
# edge.num.opt: a vector of numbers of edges for each estimated graph
# edge.opt: a list of edges for each estimated graph
# adj.mat.opt: a list of adjacency matrices for each estimated graph

loggle.cv.select.h <- function(cv.result, select.type = "all_flexible", cv.vote.thres = 0.8) {
  
  cv.score <- cv.result$cv.score
  cv.score[is.na(cv.score)] <- Inf
  cv.score[is.nan(cv.score)] <- Inf
  cv.score[is.infinite(cv.score)] <- Inf
  cv.result.fold <- cv.result$cv.result.fold
  
  L <- dim(cv.score)[1]
  D <- dim(cv.score)[2]
  K <- dim(cv.score)[3]
  cv.fold <- dim(cv.score)[4]
  p <- dim(cv.result.fold[[1]]$Omega[[L, D, K]])[1]
  
  lambda.list <- as.numeric(rownames(cv.score))
  d.list <- as.numeric(colnames(cv.score))
  
  adj.mat <- array(NA, c(p, p, K, cv.fold))
  edge.num.opt <- rep(NA, K)
  edge.opt <- vector("list", K)
  adj.mat.opt <- vector("list", K)
  
  cv.score.fold <- rowMeans(cv.score, dims = 3)
  
  d.index <- rep(0, K)
  lambda.index <- rep(0, K)
  
  if(select.type == "all_flexible") {
    
    for(k in 1:K) {
      index <- which(cv.score.fold[, , k, drop = FALSE] == min(cv.score.fold[, , k]), arr.ind = TRUE)
      index <- index[nrow(index), ]
      d.index[k] <- index[2]
      lambda.index[k] <- index[1]
    }
    
  } else if(select.type == "d_fixed") {
    
    d.index <- rep(which.min(sapply(1:D, function(d) mean(apply(cv.score.fold[, d, ], 2, min)))), K)
    lambda.index <- sapply(1:K, function(k) which.min(cv.score.fold[, d.index[1], k]))
    
  } else if(select.type == "all_fixed") {
    
    cv.score.fold.avg <- rowMeans(cv.score.fold, dims = 2)
    index <- which(cv.score.fold.avg == min(cv.score.fold.avg), arr.ind = TRUE)
    index <- index[nrow(index), ]
    d.index <- rep(index[2], K)
    lambda.index <- rep(index[1], K)
  }
  
  d.opt <- d.list[d.index]
  lambda.opt <- lambda.list[lambda.index]
  
  cv.temp <- sapply(1:cv.fold, function(i) mean(sapply(1:K, function(k) cv.score[lambda.index[k], d.index[k], k, i])))
  cv.score.opt <- mean(cv.temp)
  
  for(i in 1:cv.fold) {
    
    Omega <- cv.result.fold[[i]]$Omega
    for(k in 1:K) {
      adj.mat[, , k, i] <- as.matrix(Omega[[lambda.index[k], d.index[k], k]])
    }
  }
  
  adj.mat <- rowSums(adj.mat != 0, dims = 3) >= cv.fold * cv.vote.thres
  
  for(k in 1:K) {
    
    edge <- which(adj.mat[, , k] != 0, arr.ind = TRUE)
    edge.opt[[k]] <- edge[(edge[, 1] - edge[, 2]) > 0, , drop = FALSE]
    edge.num.opt[k] <- nrow(edge.opt[[k]])
    adj.mat.opt[[k]] <- Matrix(as.numeric(adj.mat[, , k]), p, p, sparse = TRUE)
  }
  
  result <- list(d.opt = d.opt, lambda.opt = lambda.opt, cv.score.opt = cv.score.opt, edge.num.opt = edge.num.opt, 
                 edge.opt = edge.opt, adj.mat.opt = adj.mat.opt)
  return(result)
}


# Selection function for cross validation result #######################################################################
########################################################################################################################

# Input ###
# cv.result: results from loggle.cv
# select.type: "all_flexible":optimal d and lambda can vary across time points specified by pos, 
#              "d_fixed": optimal d is fixed and optimal lambda can vary across time points specified by pos, 
#              "all_fixed": optimal d and lambda are fixed across time points specified by pos
# cv.vote.thres: an edge is kept after cv.vote if and only if it exists in no less than cv.vote.thres*cv.fold cv folds

# Output ###
# h.opt: optimal value of h
# d.opt: a vector of optimal values of d for each estimated graph
# lambda.opt: a vector of optimal values of lambda for each estimated graph
# cv.score.opt: optimal cv score (averaged over time points and cv folds)
# edge.num.opt: a vector of numbers of edges for each estimated graph
# edge.opt: a list of edges for each estimated graph
# adj.mat.opt: a list of adjacency matrices for each estimated graph

loggle.cv.select <- function(cv.result, select.type = "all_flexible", cv.vote.thres = 0.8) {
  
  H <- length(cv.result$cv.result.h)
  h.list <- sapply(1:H, function(h) cv.result$cv.result.h[[h]]$h)
  cv.select.result.h <- vector("list", H)
  cv.score.opt.h <- rep(NA, H)
  
  for(h in 1:H) {
    cv.select.result.h[[h]] <- loggle.cv.select.h(cv.result$cv.result.h[[h]], select.type, cv.vote.thres)
    cv.score.opt.h[h] <- cv.select.result.h[[h]]$cv.score.opt
  }
  
  h.opt <- h.list[which.min(cv.score.opt.h)]
  result <- cv.select.result.h[[which.min(cv.score.opt.h)]]
  result$h.opt <- h.opt

  return(result)
}


# Model refitting function given graph structures ######################################################################
########################################################################################################################

# Input ###
# X: a p by N data matrix containing observations on a time grid ranging from 0 to 1
# pos: indices of time points where graphs are estimated
# adj.mat: adjacency matrices at time points specified by pos
# h: bandwidth in kernel smoothed sample covariance/correlation matrix
# print.detail: if TRUE, print details in model refitting procedure

# Output ###
# Omega: a list of estimated precision matrices at time points specified by pos

loggle.refit <- function(X, pos, adj.mat, h = 0.8*ncol(X)^(-1/5), print.detail = TRUE) {
  
  p <- dim(X)[1]
  N <- dim(X)[2]
  K <- length(pos)
  
  if(any(!pos %in% 1:N)) {
    stop("pos must be a subset of 1, 2, ..., N!")
  }
  if(length(adj.mat) != K) {
    stop("adj.mat must have the same length as pos!")
  }
  
  cat("Generating sample covariance matrices...\n")
  Sigma <- makeCorr(X, 1:N, h, fit.corr = FALSE)$Corr
  
  cat("Estimating graphs...\n")
  
  Omega <- vector("list", K)
  
  for(k in 1:K) {
    
    edge.zero <- which(as.matrix(adj.mat[[k]]) == 0, arr.ind = TRUE)
    edge.zero <- edge.zero[(edge.zero[, 1] - edge.zero[, 2]) > 0, , drop = FALSE]
    if(nrow(edge.zero) == 0) {
      edge.zero = NULL
    }
    
    Omega[[k]] <- glasso(s = Sigma[, , pos[k]], rho = 1e-10, zero = edge.zero)$wi
    if(det(Omega[[k]]) < 0) {
      Omega[[k]] <- glasso(s = Sigma[, , pos[k]], rho = 1e-10, zero = edge.zero, thr = 5*1e-5)$wi
    }
    Omega[[k]] <- Matrix(Omega[[k]], sparse = TRUE)
    
    if(print.detail) {
      cat("Complete: t =", round((pos[k]-1) / (N-1), 2), "\n")
    }
  }
  
  return(Omega)
}


# Graph estimation function for loggle using cv.vote ###################################################################
########################################################################################################################

# Input ###
# X: a p by N data matrix containing observations on a time grid ranging from 0 to 1
# pos: indices of time points where graphs are estimated
# h: bandwidth in kernel smoothed sample covariance/correlation matrix
# d: width of neighborhood centered at each time point specified by pos
# lambda: tuning parameter of lasso penalty at each time point specified by pos
# cv.fold: number of cross-validation folds
# fit.type: likelihood: likelihood estimation, 
#           pseudo: pseudo likelihood estimation, 
#           space: sparse partial correlation estimation
# refit: if TRUE, conduct model refitting given learned graph structures
# cv.vote.thres: an edge is kept after cv.vote if and only if it exists in no less than cv.vote.thres*cv.fold cv folds
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion
# max.step: maximum steps in ADMM iteration
# detrend: if TRUE, subtract kernel weighted moving average for each variable in data matrix,
#          if FALSE, subtract overall average for each variable in data matrix
# fit.corr: if TRUE, use sample correlation matrix in model fitting,
#           if FALSE, use sample covariance matrix in model fitting
# num.thread: number of threads used in parallel computing
# print.detail: if TRUE, print details in model fitting procedure

# Output ###
# result.fold: a list of model fitting results from loggle for each cv fold
# Omega: a list of estimated precision matrices at time points specified by pos
# edge.num: a vector of numbers of edges at time points specified by pos
# edge: a list of edges at time points specified by pos

loggle.cv.vote <- function(X, pos = 1:ncol(X), h = 0.8*ncol(X)^(-1/5), d = 0.2, lambda = 0.25, cv.fold = 5, 
                           fit.type = "pseudo", refit = TRUE, cv.vote.thres = 0.8, epi.abs = 1e-5, epi.rel = 1e-3, 
                           max.step = 500, detrend = TRUE, fit.corr = TRUE, num.thread = 1, print.detail = TRUE) {
  
  p <- dim(X)[1]
  N <- dim(X)[2]
  K <- length(pos)
  
  result.fold <- vector("list", cv.fold)
  adj.mat <- array(NA, c(p, p, K, cv.fold))
  edge.num.list <- rep(NA, K)
  edge.list <- vector("list", K)
  Omega.list <- vector("list", K)
  
  for(i in 1:cv.fold) {
    
    cat("\nRunning fold", i, "out of", cv.fold, "folds...\n")
    
    pos.test <- seq(i, N, cv.fold)
    pos.train <- (1:N)[-pos.test]
    
    result.fold[[i]] <- loggle(X, pos, h, d, lambda, fit.type, refit = TRUE, epi.abs, epi.rel, max.step, detrend, 
                               fit.corr, pos.train, num.thread, print.detail)
  }
  
  cat(sprintf("\nSelecting models based on %d-fold cross-validation results...\n", cv.fold))
  
  for(i in 1:cv.fold) {
    
    Omega <- result.fold[[i]]$Omega
    for(k in 1:K) {
      adj.mat[, , k, i] <- as.matrix(Omega[[k]])
    }
  }
  
  adj.mat <- rowSums(adj.mat != 0, dims = 3) >= cv.fold * cv.vote.thres
  
  for(k in 1:K) {
    
    edge <- which(adj.mat[, , k] != 0, arr.ind = TRUE)
    edge.list[[k]] <- edge[(edge[, 1] - edge[, 2]) > 0, , drop = FALSE]
    edge.num.list[k] <- nrow(edge.list[[k]])
    Omega.list[[k]] <- Matrix(as.numeric(adj.mat[, , k]), p, p, sparse = TRUE)
  }
  
  if(refit) {
    Omega.list <- loggle.refit(X, pos, Omega.list, h)
  }
  
  result <- list(result.fold = result.fold, Omega = Omega.list, edge.num = edge.num.list, edge = edge.list)
  return(result)
}


# Cross validation function for local loggle ###########################################################################
########################################################################################################################

# Input ###
# pos: index of time point where graph is estimated
# Corr: list of kernel smoothed sample covariance/correlation matrices
# sd.X: if fit.corr = TRUE: list of standard deviations of variables
#       if fit.corr = FALSE: list of 1's
# d.list: a grid of widths of neighborhood centered at each time point specified by pos
# lambda.list: a grid of tuning parameters of lasso penalty at each time point specified by pos
# fit.type: 0: likelihood estimation, 
#           1: pseudo likelihood estimation, 
#           2: sparse partial correlation estimation
# early.stop.thres: grid search stops when the ratio between edge number and variable number exceeds early.stop.thres
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion
# max.step: maximum steps in ADMM iteration
# print.detail: if TRUE, print details in model fitting procedure

# Output ###
# Omega.list: a list of estimated precision matrices via model refitting for each combination of d and lambda
# edge.num.list: a matrix of numbers of edges for each combination of d and lambda
# edge.list: a list of edges for each combination of d and lambda

loggle.local.cv <- function(pos, Corr, sd.X, d.list, lambda.list, fit.type, early.stop.thres, epi.abs, epi.rel, 
                            max.step, print.detail) {
  
  p <- dim(Corr)[1]
  N <- dim(Corr)[3]
  D <- length(d.list)
  L <- length(lambda.list)
  
  Omega.list <- matrix(vector("list", 1), L, D)
  edge.num.list <- matrix(0, L, D)
  edge.list <- matrix(vector("list", 1), L, D)
  
  for(j in 1:D) {
    
    d <- d.list[j]
    epi.abs.d <- epi.abs[j]
    epi.rel.d <- epi.rel[j]
    
    Nd.index <- max(1, ceiling(((pos-1)/(N-1)-d)*(N-1)-1e-5)+1) : min(N, floor(((pos-1)/(N-1)+d)*(N-1)+1e-5)+1)
    Nd <- length(Nd.index)
    Nd.pos <- which(Nd.index == pos)
    Nd.pos.c <- Nd.pos - 1
    Nd.pos.l <- 1
    
    Corr.sq <- rowSums(Corr[, , Nd.index, drop = FALSE]^2, dims = 2)
    
    Z.vec <- rep(0, p*p*L)
    
    lambda <- sqrt(Nd) * lambda.list
    rho <- lambda
    
    # detect block diagonal structure
    member.index.list <- rep(0, p*L)
    no.list <- rep(0, L)
    csize.index.list <- c()
    
    for(l in L:1) {
      
      adj.mat <- (Corr.sq > lambda[l]^2)
      diag(adj.mat) <- 1
      graph <- graph.adjacency(adj.mat)
      cluster <- clusters(graph)
      member <- cluster$membership
      csize <- cluster$csize
      no <- cluster$no
      member.index <- sort(member, index.return = TRUE)$ix - 1
      csize.index <- c(0, cumsum(csize))
      
      member.index.list[(p*(l-1)+1) : (p*l)] <- member.index
      no.list[l] <- no
      csize.index.list <- c(csize.index.list, csize.index)
    }
      
    result <- .C("ADMM_lambda",
                 as.double(Corr[, , Nd.index]),
                 Z.vec = as.double(Z.vec),
                 as.integer(p),
                 as.integer(Nd),
                 as.integer(Nd.pos.c),
                 as.integer(Nd.pos.l),
                 as.integer(member.index.list),
                 as.integer(csize.index.list),
                 as.integer(no.list),
                 as.double(lambda),
                 as.integer(L),
                 as.integer(fit.type),
                 as.double(early.stop.thres),
                 as.double(rho),
                 as.double(epi.abs.d),
                 as.double(epi.rel.d),
                 as.integer(max.step)
    )
      
    Z.vec <- result$Z.vec
      
    for(l in L:1) {
      
      if(Z.vec[p*p*l] != -1) {
      
        Omega <- matrix(Z.vec[(p*p*(l-1)+1) : (p*p*l)], p, p)
        
        edge <- which(Omega != 0, arr.ind = TRUE)
        edge <- edge[(edge[, 1] - edge[, 2]) > 0, , drop = FALSE]
        edge.num <- nrow(edge)
          
        edge.zero <- which(Omega == 0, arr.ind = TRUE)
        edge.zero <- edge.zero[(edge.zero[, 1] - edge.zero[, 2]) > 0, , drop = FALSE]
        if(nrow(edge.zero) == 0) {
          edge.zero = NULL
        }
          
        Sigma <- diag(sd.X) %*% Corr[, , pos] %*% diag(sd.X)
        Omega <- glasso(s = Sigma, rho = 1e-10, zero = edge.zero)$wi
        if(det(Omega) < 0) {
          Omega <- glasso(s = Sigma, rho = 1e-10, zero = edge.zero, thr = 5*1e-5)$wi
        }
      
        Omega.list[[l, j]] <- Matrix(Omega, sparse = TRUE)
        edge.num.list[l, j] <- edge.num
        edge.list[[l, j]] <- edge
      } else {
        
        Omega.list[[l, j]] <- NA
        edge.num.list[l, j] <- NA
        edge.list[[l, j]] <- NA
      }
    }
    
    if(print.detail) {
      cat(sprintf("Complete: d = %.3f, t = %.2f\n", d, round((pos-1)/(N-1), 2)))
    }
  }
  
  result <- list(Omega.list = Omega.list, edge.num.list = edge.num.list, edge.list = edge.list)
  return(result)
}


# Cross validation function for global loggle ##########################################################################
########################################################################################################################

# Input ###
# pos: indices of time points where graphs are estimated
# Corr: list of kernel smoothed sample covariance/correlation matrices
# sd.X: if fit.corr = TRUE: list of standard deviations of variables
#       if fit.corr = FALSE: list of 1's
# lambda.list: a grid of tuning parameters of lasso penalty at each time point specified by pos
# fit.type: 0: likelihood estimation, 
#           1: pseudo likelihood estimation, 
#           2: sparse partial correlation estimation
# early.stop.thres: grid search stops when the ratio between edge number and variable number exceeds early.stop.thres
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion
# max.step: maximum steps in ADMM iteration
# print.detail: if TRUE, print details in model fitting procedure

# Output ###
# Omega.list: a list of estimated precision matrices via model refitting for each combination of lambda and time point
# edge.num.list: a matrix of numbers of edges for each combination of lambda and time point
# edge.list: a list of edges for each combination of lambda and time point

loggle.global.cv <- function(pos, Corr, sd.X, lambda.list, fit.type, early.stop.thres, epi.abs, epi.rel, max.step, 
                             print.detail) {
  
  p <- dim(Corr)[1]
  N <- dim(Corr)[3]
  K <- length(pos)
  L <- length(lambda.list)
  
  Omega.list <- array(vector("list", 1), c(L, 1, K))
  edge.num.list <- array(0, c(L, 1, K))
  edge.list <- array(vector("list", 1), c(L, 1, K))
  
  N.index.c <- 0:(N-1)
  pos.c <- pos - 1
  
  Corr.sq <- rowSums(Corr^2, dims = 2)
  
  Z.vec <- rep(0, p*p*K*L)
  
  lambda <- sqrt(N) * lambda.list
  rho <- lambda
  
  # detect block diagonal structure
  member.index.list <- rep(0, p*L)
  no.list <- rep(0, L)
  csize.index.list <- c()
  
  for(l in L:1) {
    
    adj.mat <- (Corr.sq > lambda[l]^2)
    diag(adj.mat) <- 1
    graph <- graph.adjacency(adj.mat)
    cluster <- clusters(graph)
    member <- cluster$membership
    csize <- cluster$csize
    no <- cluster$no
    member.index <- sort(member, index.return = TRUE)$ix - 1
    csize.index <- c(0, cumsum(csize))
    
    member.index.list[(p*(l-1)+1) : (p*l)] <- member.index
    no.list[l] <- no
    csize.index.list <- c(csize.index.list, csize.index)
  }
  
  result <- .C("ADMM_lambda",
               as.double(Corr),
               Z.vec = as.double(Z.vec),
               as.integer(p),
               as.integer(N),
               as.integer(pos.c),
               as.integer(K),
               as.integer(member.index.list),
               as.integer(csize.index.list),
               as.integer(no.list),
               as.double(lambda),
               as.integer(L),
               as.integer(fit.type),
               as.double(early.stop.thres),
               as.double(rho),
               as.double(epi.abs),
               as.double(epi.rel),
               as.integer(max.step)
  )
  
  Z.vec <- result$Z.vec
  
  for(l in L:1) {
    
    if(Z.vec[p*p*K*l] != -1) {
    
      Omega <- array(Z.vec[(p*p*K*(l-1)+1) : (p*p*K*l)], c(p, p, K))
      
      edge <- which(Omega[, , 1] != 0, arr.ind = TRUE)
      edge <- edge[(edge[, 1] - edge[, 2]) > 0, , drop = FALSE]
      edge.num <- nrow(edge)
        
      edge.zero <- which(Omega[, , 1] == 0, arr.ind = TRUE)
      edge.zero <- edge.zero[(edge.zero[, 1] - edge.zero[, 2]) > 0, , drop = FALSE]
      if(nrow(edge.zero) == 0) {
        edge.zero = NULL
      }
      
      Omega <- array(0, c(p, p, K))
      for(k in 1:K) {
        Sigma <- diag(sd.X) %*% Corr[, , pos[k]] %*% diag(sd.X)
        Omega[, , k] <- glasso(s = Sigma, rho = 1e-10, zero = edge.zero)$wi
        if(det(Omega[, , k]) < 0) {
          Omega[, , k] <- glasso(s = Sigma, rho = 1e-10, zero = edge.zero, thr = 5*1e-5)$wi
        }
      }
      
      edge.num.list[l, 1, ] <- edge.num
    
      for(k in 1:K) {
      
        Omega.list[[l, 1, k]] <- Matrix(Omega[, , k], sparse = TRUE)
        edge.list[[l, 1, k]] <- edge
      }
    } else {
      
      edge.num.list[l, 1, ] <- NA
      for(k in 1:K) {
        Omega.list[[l, 1, k]] <- NA
        edge.list[[l, 1, k]] <- NA
      }
    }
  }
  
  if(print.detail) {
    cat("Complete: d = 1", "\n")
  }
  
  result <- list(Omega.list = Omega.list, edge.num.list = edge.num.list, edge.list = edge.list)
  return(result)
}


# Cross validation function for local&global loggle ####################################################################
########################################################################################################################

# Input ###
# X: a p by N data matrix containing observations on a time grid ranging from 0 to 1
# pos.train: indices of time points used to generate sample covariance/correlation matrix
# pos: indices of time points where graphs are estimated
# h:  bandwidth in kernel smoothed sample covariance/correlation matrix
# d.list: a grid of widths of neighborhood centered at each time point specified by pos
# lambda.list: a grid of tuning parameters of lasso penalty at each time point specified by pos
# fit.type: 0: likelihood estimation, 
#           1: pseudo likelihood estimation, 
#           2: sparse partial correlation estimation
# early.stop.thres: grid search stops when the ratio between edge number and variable number exceeds early.stop.thres
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion
# max.step: maximum steps in ADMM iteration
# fit.corr: if TRUE, use sample correlation matrix in model fitting,
#           if FALSE, use sample covariance matrix in model fitting
# num.thread: number of threads used in parallel computing
# print.detail: if TRUE, print details in model fitting procedure

# Output ###
# Omega: a list of estimated precision matrices via model refitting for each combination of d, lambda and time point
# edge.num: an array of numbers of edges for each combination of d, lambda and time point
# edge: a list of edges for each combination of d, lambda and time point

loggle.combine.cv <- function(X, pos.train, pos, h, d.list, lambda.list, fit.type, early.stop.thres, epi.abs, epi.rel, 
                              max.step, fit.corr, num.thread, print.detail) {
  
  p <- dim(X)[1]
  N <- dim(X)[2]
  K <- length(pos)
  D <- length(d.list)
  L <- length(lambda.list)
  
  if(fit.corr) {
    cat("Generating sample correlation matrices for training dataset...\n")
  } else {
    cat("Generating sample covariance matrices for training dataset...\n")
  }
  result.Corr <- makeCorr(X, pos.train, h, fit.corr)
  Corr <- result.Corr$Corr
  sd.X <- result.Corr$sd.X
  rm(result.Corr)
  
  cat("Estimating graphs for training dataset...\n")
  
  if(d.list[1] == 1) {

    result <- loggle.global.cv(pos, Corr, sd.X, lambda.list, fit.type, early.stop.thres, epi.abs, epi.rel, max.step, 
                               print.detail)
    result <- list(Omega = result$Omega.list, edge.num = result$edge.num.list, edge = result$edge.list)
    
  } else {
    
    Omega.list <- array(vector("list", 1), c(L, D, K))
    edge.num.list <- array(NA, c(L, D, K))
    edge.list <- array(vector("list", 1), c(L, D, K))
    
    if(d.list[D] == 1) {
      
      if(num.thread > 1) {
        
        result <- foreach(k=1:K, .combine="list", .multicombine=TRUE, .maxcombine=K, .export=c("loggle.local.cv")) %dopar%
          loggle.local.cv(pos[k], Corr, sd.X, d.list[-D], lambda.list, fit.type, early.stop.thres, epi.abs[-D], 
                          epi.rel[-D], max.step, print.detail)
      } else {
        
        result <- vector("list", K)
        for(k in 1:K) {
          result[[k]] <- loggle.local.cv(pos[k], Corr, sd.X, d.list[-D], lambda.list, fit.type, early.stop.thres, 
                                         epi.abs[-D], epi.rel[-D], max.step, print.detail)
        }
      }
      
      for(k in 1:K) {
        
        result.k <- result[[k]]
        Omega.list[, -D, k] <- result.k$Omega.list
        edge.num.list[, -D, k] <- result.k$edge.num.list
        edge.list[, -D, k] <- result.k$edge.list
      }
      
      result <- loggle.global.cv(pos, Corr, sd.X, lambda.list, fit.type, early.stop.thres, epi.abs[D], epi.rel[D], 
                                 max.step, print.detail)
      
      Omega.list[, D, ] <- result$Omega.list
      edge.num.list[, D, ] <- result$edge.num.list
      edge.list[, D, ] <- result$edge.list
      
    } else {
      
      if(num.thread > 1) {
        
        result <- foreach(k=1:K, .combine="list", .multicombine=TRUE, .maxcombine=K, .export=c("loggle.local.cv")) %dopar%
          loggle.local.cv(pos[k], Corr, sd.X, d.list, lambda.list, fit.type, early.stop.thres, epi.abs, epi.rel, 
                          max.step, print.detail)
      } else {
        
        result <- vector("list", K)
        for(k in 1:K) {
          result[[k]] <- loggle.local.cv(pos[k], Corr, sd.X, d.list, lambda.list, fit.type, early.stop.thres, epi.abs, 
                                         epi.rel, max.step, print.detail)
        }
      }
      
      for(k in 1:K) {
        
        result.k <- result[[k]]
        Omega.list[, , k] <- result.k$Omega.list
        edge.num.list[, , k] <- result.k$edge.num.list
        edge.list[, , k] <- result.k$edge.list
      }
    }
    
    result <- list(Omega = Omega.list, edge.num = edge.num.list, edge = edge.list)
  }

  return(result)  
}