
# Covariance/Correlation matrix generation function ###########################################################################
###############################################################################################################################

# Input ###
# X: a p by N matrix containing list of observations
# pos: position of observations used to generate correlation matrices
# h: bandwidth in kernel function used to generate correlation matrices
# fit.corr: whether to use sample correlation matrix rather than sample covariance matrix in model fitting

# Output ###
# Corr: list of kernel estimators of correlation matrices
# sd.X: list of standard deviations of variables

makeCorr <- function(X, pos, h, fit.corr) {
  
  p <- nrow(X)
  N <- ncol(X)
  L <- length(pos)
  
  sd.X <- rep(NA, p)
  
  if(fit.corr) {
    for(i in 1:p) {
      sd.X[i] <- sd(X[i, ])
      X[i, ] <- scale(X[i, ])
    }
  } else {
    for(i in 1:p) {
      sd.X[i] <- 1
      X[i, ] <- scale(X[i, ], scale = FALSE)
    }
  }
  
  Corr <- array(0, c(p, p, N))
  
  for(i in 1:N) {
    Kh <- pmax(3/4 * (1 - ((pos - i) / ((N - 1) * h)) ^ 2), 0)
    omega <- Kh / sum(Kh)
    index <- which(omega != 0)
    Corr[, , i] <- X[, pos[index]] %*% diag(omega[index]) %*% t(X[, pos[index]])
  }
  
  result <- new.env()
  result$Corr <- Corr
  result$sd.X <- sd.X
  result <- as.list(result)
  
  return(result)
}


# Cross validation function for local LGGM ####################################################################################
###############################################################################################################################

# Input ###
# pos: position of time point where graph is estimated
# Corr: list of kernel estimators of correlation matrices
# sd.X: list of standard deviations of variables
# fit.type: 0: graphical Lasso estimation, 1: pseudo likelihood estimation, 3: sparse partial correlation estimation
# refit.type: 0: likelihood estimation, 1: pseudo likelihood estimation
# d.list: list of widths of neighborhood
# lambda.list: list of tuning parameters of Lasso penalty
# cv.thres: grid search stops when number of detected edges exceeds cv.thres times number of nodes
# epi.abs: list of absolute tolerances in ADMM stopping criterion
# epi.rel: list of relative tolerances in ADMM stopping criterion

# Output ###
# Omega: L (number of lambda's) by D (number of d's) list of estimated precision matrices
# Omega.rf: L by D list of refitted precision matrices
# edge.num: L by D list of edge numbers
# edge: L by D list of edges

LGGM.local.cv <- function(pos, Corr, sd.X, fit.type, refit.type, d.list, lambda.list, cv.thres, epi.abs, epi.rel) {
  
  p <- dim(Corr)[1]
  N <- dim(Corr)[3]
  D <- length(d.list)
  L <- length(lambda.list)
  
  Omega.list <- matrix(vector("list", 1), L, D)
  Omega.rf.list <- matrix(vector("list", 1), L, D)
  edge.num.list <- matrix(0, L, D)
  edge.list <- matrix(vector("list", 1), L, D)
  
  for(j in 1:D) {
    
    d <- d.list[j]
    epi.abs.d <- epi.abs[j]
    epi.rel.d <- epi.rel[j]
    
    Nd.index <- max(1, ceiling(((pos-1)/(N-1) - d) * (N-1) - 1e-5) + 1) : min(N, floor(((pos-1)/(N-1) + d) * (N-1) + 1e-5) + 1)
    Nd <- length(Nd.index)
    Nd.pos <- which(Nd.index == pos)
    Nd.pos.c <- Nd.pos - 1
    Nd.pos.l <- 1
    
    Corr.sq <- apply(Corr[, , Nd.index] ^ 2, c(1, 2), sum)
    
    Z.vec <- rep(0, p*p*L)
    Z.pos.vec <- rep(0, p*p*L)
    
    lambda <- sqrt(Nd) * lambda.list
    rho <- lambda
    
    # detect block diagonal structure
    member.index.list <- rep(0, p*L)
    no.list <- rep(0, L)
    csize.index.list <- c()
    
    for(l in L:1) {
      
      adj.mat <- (Corr.sq > lambda[l] ^ 2)
      diag(adj.mat) <- 1
      graph <- graph.adjacency(adj.mat)
      cluster <- clusters(graph)
      member <- cluster$membership
      csize <- cluster$csize
      no <- cluster$no
      member.index <- sort(member, index.return=T)$ix - 1
      csize.index <- c(0, cumsum(csize))
      
      member.index.list[(p*(l-1)+1) : (p*l)] <- member.index
      no.list[l] <- no
      csize.index.list <- c(csize.index.list, csize.index)
    }
      
    result <- .C("ADMM_lambda", 
                  as.integer(p),
                  as.integer(member.index.list),
                  as.integer(csize.index.list),
                  as.integer(no.list),
                  as.integer(Nd),
                  as.integer(Nd.pos.c),
                  as.integer(Nd.pos.l),
                  as.double(Corr[, , Nd.index]),
                  as.double(sd.X),
                  Z.ini.vec = as.double(Z.vec),
                  Z.pos.vec = as.double(Z.pos.vec),
                  as.integer(L),
                  as.double(lambda),
                  as.double(rho),
                  as.double(epi.abs.d),
                  as.double(epi.rel.d),
                  as.integer(pseudo.fit),
                  as.integer(pseudo.refit),
                  as.double(cv.thres)
    )
      
    Z.vec <- result$Z.vec
    Z.pos.vec <- result$Z.pos.vec
      
    for(l in L:1){
      
      Omega <- matrix(Z.vec[(p*p*(l-1)+1) : (p*p*l)], p, p)
      Omega.rf <- matrix(Z.pos.vec[(p*p*(l-1)+1) : (p*p*l)], p, p)
        
      edge <- which(Omega.rf != 0, arr.ind = T)
      edge <- edge[(edge[, 1] - edge[, 2]) > 0, , drop = F]
      edge.num <- nrow(edge)
        
      Omega.list[[l, j]] <- Matrix(Omega, sparse = T)
      Omega.rf.list[[l, j]] <- Matrix(Omega.rf, sparse = T)
      edge.num.list[l, j] <- edge.num
      edge.list[[l, j]] <- edge
    }
    
    cat("Complete: d =", d, "t =", round((pos-1)/(N-1), 2), "\n")
    
    result <- new.env()
    result$Omega.list <- Omega.list
    result$Omega.rf.list <- Omega.rf.list
    result$edge.num.list <- edge.num.list
    result$edge.list <- edge.list
    result <- as.list(result)
  }
  
  return(result)
}


# Cross validation function for global LGGM ###################################################################################
###############################################################################################################################

# Input ###
# pos: position of time points where graphs are estimated
# Corr: list of kernel estimators of correlation matrices
# sd.X: list of standard deviations of variables
# fit.type: 0: graphical Lasso estimation, 1: pseudo likelihood estimation, 3: sparse partial correlation estimation
# refit.type: 0: likelihood estimation, 1: pseudo likelihood estimation
# lambda.list: list of tuning parameters of Lasso penalty
# cv.thres: grid search stops when number of detected edges exceeds cv.thres times number of nodes
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion

# Output ###
# Omega: L (number of lambda's) by K (number of time points) list of estimated precision matrices
# Omega.rf: L by K list of refitted precision matrices
# edge.num: L by K list of edge numbers
# edge: L by K list of edges

LGGM.global.cv <- function(pos, Corr, sd.X, fit.type, refit.type, lambda.list, cv.thres, epi.abs, epi.rel) {
  
  p <- dim(Corr)[1]
  N <- dim(Corr)[3]
  K <- length(pos)
  L <- length(lambda.list)
  
  Omega.list <- array(vector("list", 1), L, 1, K)
  Omega.rf.list <- array(vector("list", 1), L, 1, K)
  edge.num.list <- array(0, c(L, 1, K))
  edge.list <- array(vector("list", 1), L, 1, K)
  
  N.index.c <- 0 : (N-1)
  pos.c <- pos - 1
  
  Corr.sq <- apply(Corr ^ 2, c(1, 2), sum)
  
  Z.vec <- rep(0, p*p*K*L)
  Z.pos.vec <- rep(0, p*p*K*L)
  
  lambda <- sqrt(N) * lambda.list
  rho <- lambda
  
  # detect block diagonal structure
  member.index.list <- rep(0, p*L)
  no.list <- rep(0, L)
  csize.index.list <- c()
  
  for(l in L:1) {
    
    adj.mat <- (Corr.sq > lambda[l] ^ 2)
    diag(adj.mat) <- 1
    graph <- graph.adjacency(adj.mat)
    cluster <- clusters(graph)
    member <- cluster$membership
    csize <- cluster$csize
    no <- cluster$no
    member.index <- sort(member, index.return = T)$ix - 1
    csize.index <- c(0, cumsum(csize))
    
    member.index.list[(p*(l-1)+1) : (p*l)] <- member.index
    no.list[l] <- no
    csize.index.list <- c(csize.index.list, csize.index)
  }
  
  result <- .C("ADMM_lambda", 
              as.integer(p),
              as.integer(member.index.list),
              as.integer(csize.index.list),
              as.integer(no.list),
              as.integer(N),
              as.integer(pos.c),
              as.integer(K),
              as.double(Corr),
              as.double(sd.X),
              Z.ini.vec = as.double(Z.vec),
              Z.pos.vec = as.double(Z.pos.vec),
              as.integer(L),
              as.double(lambda),
              as.double(rho),
              as.double(epi.abs),
              as.double(epi.rel),
              as.integer(pseudo.fit),
              as.integer(pseudo.refit),
              as.double(cv.thres)
  )
  
  Z.vec <- result$Z.vec
  Z.pos.vec <- result$Z.pos.vec
  
  for(l in L:1) {
    
    Omega <- array(Z.vec[(p*p*K*(l-1)+1) : (p*p*K*l)], c(p, p, K))
    Omega.rf <- array(Z.pos.vec[(p*p*K*(l-1)+1) : (p*p*K*l)], c(p, p, K))
    
    edge <- which(Omega.rf[, , 1] != 0, arr.ind = T)
    edge <- edge[(edge[, 1] - edge[, 2]) > 0, , drop = F]
    edge.num <- nrow(edge)
  
    edge.num.list[l, 1, ] <- edge.num
    
    for(k in 1:K) {
      
      Omega.list[[l, 1, k]] <- Matrix(Omega[, , k], sparse = T)
      Omega.rf.list[[l, 1, k]] <- Matrix(Omega.rf[, , k], sparse = T)
      edge.list[[l, 1, k]] <- edge
    }
  }
  
  cat("Complete: d = 1", "\n")
  
  result <- new.env()
  result$Omega.list <- Omega.list
  result$Omega.rf.list <- Omega.rf.list
  result$edge.num.list <- edge.num.list
  result$edge.list <- edge.list
  result <- as.list(result)
  
  return(result)
}


# Cross validation function for local&global LGGM #############################################################################
###############################################################################################################################

# Input ###
# X: a p by N matrix containing list of observations
# pos: position of time points where graphs are estimated
# fit.type: 0: graphical Lasso estimation, 1: pseudo likelihood estimation, 3: sparse partial correlation estimation
# refit.type: 0: likelihood estimation, 1: pseudo likelihood estimation
# h: bandwidth in kernel function used to generate correlation matrices
# d.list: list of widths of neighborhood
# lambda.list: list of tuning parameters of Lasso penalty
# cv.thres: grid search stops when number of detected edges exceeds cv.thres times number of nodes
# epi.abs: list of absolute tolerances in ADMM stopping criterion
# epi.rel: list of relative tolerances in ADMM stopping criterion
# fit.corr: whether to use sample correlation matrix rather than sample covariance matrix in model fitting
# num.thread: number of threads

# Output ###
# Omega: L (number of lambda's) by D (number of d's) by K (number of time points) list of estimated precision matrices
# Omega.rf: L by D by K list of refitted precision matrices
# edge.num: L by D by K list of edge numbers
# edge: L by D by K list of edges

LGGM.combine.cv <- function(X, pos, fit.type, refit.type, h, d.list, lambda.list, cv.thres, epi.abs, epi.rel, fit.corr, num.thread) {
  
  p <- dim(X)[1]
  N <- dim(X)[2]
  K <- length(pos)
  D <- length(d.list)
  L <- length(lambda.list)
  
  cat("Generating sample covariance/correlation matrices...\n")
  result.Corr <- makeCorr(X, 1:N, h, fit.corr)
  Corr <- result.Corr$Corr
  sd.X <- result.Corr$sd.X
  rm(result.Corr)
  
  cat("Estimating graphs...\n")
  
  if(length(d.list) == 1 & d.list == 1) {

    result <- LGGM.global.cv(pos, Corr, sd.X, fit.type, refit.type, lambda.list, cv.thres, epi.abs, epi.rel)
    
  } else {
    
    Omega.list <- array(vector("list", 1), c(L, D, K))
    Omega.rf.list <- array(vector("list", 1), c(L, D, K))
    edge.num.list <- array(vector("list", 1), c(L, D, K))
    edge.list <- array(0, c(L, D, K))
    
    if(d.list[D] == 1) {
      
      registerDoParallel(num.thread)
      
      result <- foreach(k = 1:K, .combine = "list", .multicombine = TRUE, .maxcombine = K, .export = c("LGGM.local.cv")) %dopar%
        LGGM.local.cv(pos[k], Corr, sd.X, fit.type, refit.type, d.list[-D], lambda.list, cv.thres, epi.abs[-D], epi.rel[-D])
      
      stopImplicitCluster()
      
      for(k in 1:K) {
        
        result.k <- result[[k]]
        
        Omega.list[, -D, k] <- result.k$Omega.list
        Omega.rf.list[, -D, k] <- result.k$Omega.rf.list
        edge.num.list[, -D, k] <- result.k$edge.num.list
        edge.list[, -D, k] <- result.k$edge.list
      }
      
      result <- LGGM.global.cv(pos, Corr, Sigma, fit.type, refit.type, lambda.list, cv.thres, epi.abs[D], epi.rel[D])
      
      Omega.list[, D, ] <- result.k$Omega.list
      Omega.rf.list[, D, ] <- result.k$Omega.rf.list
      edge.num.list[, D, ] <- result.k$edge.num.list
      edge.list[, D, ] <- result.k$edge.list
      
    } else {
      
      registerDoParallel(num.thread)
      
      result <- foreach(k = 1:K, .combine = "list", .multicombine = TRUE, .maxcombine = K, .export = c("LGGM.local.cv")) %dopar%
        LGGM.local.cv(pos[k], Corr, sd.X, fit.type, refit.type, d.list, lambda.list, cv.thres, epi.abs, epi.rel)
      
      stopImplicitCluster()
      
      for(k in 1:K) {
        
        result.k <- result[[k]]
        
        Omega.list[, , k] <- result.k$Omega.list
        Omega.rf.list[, , k] <- result.k$Omega.rf.list
        edge.num.list[, , k] <- result.k$edge.num.list
        edge.list[, , k] <- result.k$edge.list
      }
      
      result <- new.env()
      result$Omega.list <- Omega.list
      result$Omega.rf.list <- Omega.rf.list
      result$edge.num.list <- edge.num.list
      result$edge.list <- edge.list
      result <- as.list(result)
    }
  }

  return(result)  
}


#LGGM.cv.select###############################################################################################################################################################################################################################################

LGGM.cv.select = function(cv.result, select.mode = "all_flexible", cv.vote.thres = 0.8){
  
  cv.score = cv.result$cv.score; cv.result.list = cv.result$cv.result.list
  
  p = dim(cv.result.list[[1]]$Omega.rf.list[[1,1,1]])[1]
  L = dim(cv.score)[1]; D = dim(cv.score)[2]; K = dim(cv.score)[3]; fold = length(cv.result.list)
  
  lambda.list = as.numeric(colnames(cv.score)); d.list = as.numeric(rownames(cv.score))
  
  Omega.edge.list.min = array(0, c(p, p, K, fold)); edge.num.list.min = rep(0, K); edge.list.min = vector("list", K)
  
  cv.score.fold = apply(cv.score, c(1,2,3), sum)
  
  lambda.index = rep(0, K); d.index = rep(0, K)
  
  if(select.mode = "all_flexible"){
    
    for(k in 1:K){
      index = which(cv.score.fold[, , k] == min(cv.score.fold[, , k]), arr.ind = T)
      if(nrow(index)>1){index = index[nrow(index), ]}
      d.index[k] = index[2]
      lambda.index[k] = index[1]
    }
  }
  if(select.mode = "d_fixed"){
    
    d.index = rep(which.min(sapply(1:D, function(d) sum(apply(cv[, d, ], 2, min)))), K)
    lambda.index = sapply(1:K, function(k) which.min(cv[, index.d, k]))
  }
  if(select.mode = "all_fixed"){
    
    cv.score.fold.avg = apply(cv.score.fold, c(1, 2), mean)
    index = which(cv.score.fold.avg == min(cv.score.fold.avg), arr.ind = T)
    d.index = rep(index[2], K)
    lambda.index = rep(index[1], K)
  }
  
  d.min = d.list[d.index]
  lambda.min = lambda.list[lambda.index]
  
  cv.temp = sapply(1:fold, function(i) mean(sapply(1:K, function(k) cv.score.fold[lambda.index[k], d.index[k], k, i])))
  cv.score.min = mean(cv.temp)
  cv.score.min.sd = sd(cv.temp)/sqrt(fold)
  
  for(i in 1:fold){
    
    Omega.rf.list = cv.result.list[[i]]$Omega.rf.list
    for(k in 1:K){Omega.edge.list.min[, , k, i] = as.matrix(Omega.rf.list[[lambda.index[k], d.index[k], k]])}
  }
  
  Omega.edge.list.min = (apply(Omega.edge.list.min!=0, c(1,2,3), sum)>=fold*cv.vote.thres)
  
  for(k in 1:K){
    
    edge = which(Omega.edge.list.min[, , k]!=0, arr.ind=T); edge.list.min[[k]] = edge[(edge[, 1] - edge[, 2])>0, , drop = F]
    edge.num.list.min[k] = nrow(edge)
  }
  
  result = new.env()
  result$d.min = d.min
  result$lambda.min = lambda.min
  result$cv.score.min = cv.score.min
  result$cv.score.min.sd = cv.score.min.sd
  result$Omega.edge.list.min = Omega.edge.list.min
  result$edge.num.list.min = edge.num.list.min
  result$edge.list.min = edge.list.min
  result = as.list(result)
  
  return(result)
}


#cross validation selection###################################################################################################################################################################################################################################

#Input###
#fold: number of folds in cross validation
#corr.ind: 0: use Sigma in model fitting, 1: use Corr in model fitting

#Output###
#result: a list of length=fold containing the results from simu.lglasso
#cv: L by D by K by fold by 2 cv scores (the last dimension corresponds to whether we use the adjusted bandwidth h or not)

LGGM.cv = function(X, pos = 1:ncol(X), fit.type = "glasso", refit.type = "likelihood", h = 0.8*ncol(X)^(-1/5), d.list, lambda.list, num.fold = 5, cv.thres = nrow(X), return.select = TRUE, select.mode = "all_flexible", cv.vote.thres = 0.8, epi.abs = 1e-5, epi.rel = 1e-3, corr = TRUE, h.correct = TRUE, num.thread = 1){
  
  p = dim(X)[1]; N = dim(X)[2]; K = length(pos); D = length(d.list); L = length(lambda.list)
  
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
  if(!length(d) %in% c(1, K)) {
    stop("d must be a scalar or a vector of the same length as 'pos'!")
  }
  if(!length(lambda) %in% c(1, K)) {
    stop("lambda must be a scalar or a vector of the same length as 'pos'!")
  }
  
  d.list.global <- d.list >= 0.5
  if(sum(d.list.global) > 0) {
    d.list.prev <- d.list
    d.list <- c(d.list[!d.list.global], 1)
    cat("Convert d.list:", d.list.prev, "to:", d.list, "\n")
  }
  d.list <- sort(d.list)
  
  lambda.list <- sort(lambda.list)
  
  
  if(length(epi.abs) == 1){
    epi.abs = rep(epi.abs, D)
  }
  if(length(epi.rel) == 1){
    epi.rel = rep(epi.rel, D)
  }
  
  if(h.correct == TRUE){
    h.test = h*(fold-1)^(1/5)
  }else{
    h.test = h
  }
  
  cv.score = array(0, c(L, D, K, fold)); rownames(cv.score) = lambda.list; colnames(cv.score) = d.list
  cv.result.list = vector("list", fold)
  
  for(i in 1:num.fold){
    
    cat("\nRunning fold", i, "out of", num.fold, "fold...\n")
    
    pos.test = seq(i, N, fold); pos.train = (1:N)[-pos.test]
    
    Sigma.test = makeCorr(X, pos.test, h.test, fit.corr = FALSE)$Corr
    
    result.i = LGGM.combine.cv(X, pos, fit.type, refit.type, h, d.list, lambda.list, cv.thres, epi.abs, epi.rel, fit.corr, num.thread)
    cv.result.list[[i]] = result.i
    
    for(d in 1:D){for(l in 1:L){for(k in 1:K){
      
      Omega.rf = as.matrix(result.i$Omega.rf.list[[l, d, k]])
      cv.score[l, d, k, i] = sum(c(t(Sigma.test[, , pos[k]]))*c(Omega.rf)) - log(det(Omega.rf))
    }}}
    
  }
  
  cv.result = new.env()
  cv.result$cv.score = cv.score
  cv.result$cv.result.list = cv.result.list
  cv.result = as.list(cv.result)
  
  if(return.select == TRUE){
    
    cv.select.result = LGGM.cv.select(cv.result, select.mode, cv.vote.thres)
    cv.result = c(cv.result, cv.select.result)
  }
  
  return(cv.result)
}


#LGGM.cv.h####################################################################################################################################################################################################################################################

LGGM.cv.h = function(X, pos.prop = 0.01, fit.type = "glasso", refit.type = "likelihood", h.list = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35), d.list = c(0, 0.01, 0.05, 0.15, 0.25, 0.35, 1), lambda.list = c(0.15, 0.2, 0.25, 0.3), fold = 5, cv.thres = 1, return.select = TRUE, select.mode = "all_flexible", cv.vote.thres = 0.8, epi.abs = 1e-4, epi.rel = 1e-2, corr = TRUE, h.correct = TRUE, num.core = 1){
  
  N = dim(X)[2]; pos = round(seq(0.02, 0.98, length=round(pos.prop*(N-1)+1))*(N-1)+1); H = length(h.list)
  
  cv.score.min.h = rep(0, H)
  
  for(h in 1:H){
    
    cv.result.h = LGGM.cv(X, pos, fit.type, refit.type, h.list[h], d.list, lambda.list, fold, cv.thres, return.select, select.mode, cv.vote.thres, epi.abs, epi.rel, corr, h.correct, num.core)
    cv.score.min.h[h] = cv.result.h$cv.score.min
  }
  
  h.min = h.list[which.min(cv.score.min.h)]
  
  return(h.min)
}