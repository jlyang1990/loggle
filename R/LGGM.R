
# Correlation matrix generation function ######################################################################################
###############################################################################################################################

# Input ###
# X: a p by N matrix containing list of observations
# pos: position of observations used to generate correlation matrices
# h: bandwidth in kernel function used to generate correlation matrices

# Output ###
# Corr: list of kernel estimators of correlation matrices
# sd.X: list of standard deviations of variables

makeCorr <- function(X, pos, h) {
  
  p <- nrow(X)
  N <- ncol(X)
  L <- length(pos)
  
  sd.X <- rep(NA, p)
  
  for(i in 1:p) {
    sd.X[i] <- sd(X[i, ])
    X[i, ] <- scale(X[i, ])
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
  result = as.list(result)
  
  return(result)
}


# Graph estimation function for local LGGM ####################################################################################
###############################################################################################################################

# Input ###
# pos: position of time point where graph is estimated
# Corr: list of kernel estimators of correlation matrices
# sd.X: list of standard deviations of variables
# fit.type: 0: graphical Lasso estimation, 1: pseudo likelihood estimation, 3: sparse partial correlation estimation
# refit.type: 0: likelihood estimation, 1: pseudo likelihood estimation
# d: width of neighborhood
# lambda: tuning parameter of Lasso penalty
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion

# Output ###
# Omega: D (number of d's) by L (number of lambda's) list of estimated precision matrices
# Omega.rf: D by L list of refitted precision matrices
# edge.num: D by L list of edge numbers
# edge: D by L list of edges

LGGM.local <- function(pos, Corr, sd.X, fit.type, refit.type, d, lambda, epi.abs, epi.rel) {
  
  p <- dim(Corr)[1]
  N <- dim(Corr)[3]
  
  Nd.index <- max(1, ceiling(((pos-1)/(N-1) - d) * (N-1) - 1e-5) + 1) : min(N, floor(((pos-1)/(N-1) + d) * (N-1) + 1e-5) + 1)
  Nd <- length(Nd.index)
  Nd.pos <- which(Nd.index == pos)
  Nd.pos.c <- Nd.pos - 1
  Nd.pos.l <- 1
    
  Corr.sq <- apply(Corr[, , Nd.index] ^ 2, c(1, 2), sum)
  
  Z.vec <- rep(0, p*p)
  Z.pos.vec <- rep(0, p*p)
  edge.num <- 0
    
  lambda <- sqrt(Nd) * lambda
  rho <- lambda
    
  #detect block diagonal structure
  adj.mat <- (Corr.sq > lambda ^ 2)
  diag(adj.mat) <- 1
  graph <- graph.adjacency(adj.mat)
  cluster <- clusters(graph)
  member <- cluster$membership
  csize <- cluster$csize
  no <- cluster$no
  member.index <- sort(member, index.return = T)$ix - 1
  csize.index <- c(0, cumsum(csize))
    
  result <- .C("ADMM_simple",
               as.integer(p),
               as.integer(member.index),
               as.integer(csize.index),
               as.integer(no),
               as.integer(Nd),
               as.integer(Nd.pos.c),
               as.integer(Nd.pos.l),
               as.double(Corr[, , Nd.index]),
               as.double(sd.X),
               Z.vec = as.double(Z.vec),
               Z.pos.vec = as.double(Z.pos.vec),
               as.double(lambda),
               as.double(rho),
               as.double(epi.abs),
               as.double(epi.rel),
               as.integer(fit.type),
               as.integer(refit.type),
               edge.num = as.integer(edge.num)
               )
  
  Omega <- Matrix(result$Z.vec, p, p, sparse = T)
  Omega.rf <- Matrix(result$Z.pos.vec, p, p, sparse = T)
  edge.num <- result$edge.num
  edge <- which(Omega.rf != 0, arr.ind = T)
  edge <- edge[(edge[, 1] - edge[, 2]) > 0, , drop = F]
    
  cat("complete: t =", round((pos-1) / (N-1), 2), "\n")
  
  result <- new.env()
  result$Omega <- Omega
  result$Omega.rf <- Omega.rf
  result$edge.num <- edge.num
  result$edge <- edge
  result <- as.list(result)
  
  return(result)
}


# Graph estimation function for global LGGM ###################################################################################
###############################################################################################################################

# Input ###
# pos: position of time points where graphs are estimated
# Corr: list of kernel estimators of correlation matrices
# sd.X: list of standard deviations of variables
# fit.type: 0: graphical Lasso estimation, 1: pseudo likelihood estimation, 3: sparse partial correlation estimation
# refit.type: 0: likelihood estimation, 1: pseudo likelihood estimation
# d: width of neighborhood
# lambda: tuning parameter of Lasso penalty
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion

# Output ###
# Omega: D (number of d's) by L (number of lambda's) by K (number of time points) list of estimated precision matrices
# Omega.rf: D by L by K list of refitted precision matrices
# edge.num: D by L by K list of edge numbers
# edge: D by L by K list of edges

LGGM.global <- function(pos, Corr, sd.X, fit.type, refit.type, lambda, epi.abs, epi.rel) {
  
  p <- dim(Corr)[1]
  N <- dim(Corr)[3]
  K <- length(pos)
  
  N.index.c <- 0:(N-1)
  pos.c <- pos - 1
  
  Corr.sq <- apply(Corr ^ 2, c(1, 2), sum)
  
  Z.vec <- rep(0, p*p*K)
  Z.pos.vec <- rep(0, p*p*K)
  edge.num <- 0
  
  lambda <- sqrt(N) * lambda
  rho <- lambda
  
  #detect block diagonal structure
  adj.mat <- (Corr.sq > lambda ^ 2)
  diag(adj.mat) <- 1
  graph <- graph.adjacency(adj.mat)
  cluster <- clusters(graph)
  member <- cluster$membership
  csize <- cluster$csize
  no <- cluster$no
  member.index <- sort(member, index.return = T)$ix - 1
  csize.index <- c(0, cumsum(csize))
  
  result <- .C("ADMM_simple",
               as.integer(p),
               as.integer(member.index),
               as.integer(csize.index),
               as.integer(no),
               as.integer(N),
               as.integer(pos.c),
               as.integer(K),
               as.double(Corr),
               as.double(sd.X),
               Z.vec = as.double(Z.vec),
               Z.pos.vec = as.double(Z.pos.vec),
               as.double(lambda),
               as.double(rho),
               as.double(epi.abs),
               as.double(epi.rel),
               as.integer(fit.type),
               as.integer(refit.type),
               edge.num = as.integer(edge.num)
  )
  
  Omega.list <- sapply(1:K, function(k) Matrix(result$Z.vec[(p*p*(k-1) + 1) : (p*p*k)], p, p, sparse = T))
  Omega.rf.list <- sapply(1:K, function(k) Matrix(result$Z.pos.vec[(p*p*(k-1) + 1) : (p*p*k)], p, p, sparse = T))
  edge.num.list <- rep(result$edge.num, K)
  edge <- which(Omega.rf.list[[1]] != 0, arr.ind = T)
  edge <- edge[(edge[, 1] - edge[, 2]) > 0, , drop = F]
  edge.list <- rep(list(edge), K)
  
  result <- new.env()
  result$Omega.list <- Omega.list
  result$Omega.rf.list <- Omega.rf.list
  result$edge.num.list <- edge.num.list
  result$edge.list <- edge.list
  result <- as.list(result)
  
  return(result)
}


# Graph estimation function for LGGM ##########################################################################################
###############################################################################################################################

# Input ###
# X: a p by N matrix containing list of observations
# pos: position of time points where graphs are estimated
# fit.type: 0: graphical Lasso estimation, 1: pseudo likelihood estimation, 3: sparse partial correlation estimation
# refit.type: 0: likelihood estimation, 1: pseudo likelihood estimation
# h: bandwidth in kernel function used to generate correlation matrices
# d: list of widths of neighborhood
# lambda: list of tuning parameters of Lasso penalty
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion
# fit.corr: whether to use sample correlation matrix rather than sample covariance matrix in model fitting
#num.thread: number of threads

# Output ###
# Omega: D (number of d's) by L (number of lambda's) by K (number of time points) list of estimated precision matrices
# Omega.rf: D by L by K list of refitted precision matrices
# edge.num: D by L by K list of edge numbers
# edge: D by L by K list of edges

LGGM <- function(X, pos = 1:ncol(X), fit.type = "glasso", refit.type = "likelihood", h = 0.8*ncol(X)^(-1/5), d, lambda, epi.abs = 1e-5, epi.rel = 1e-3, fit.corr = TRUE, num.thread = 1) {
  
  p <- dim(X)[1]
  N <- dim(X)[2]
  K <- length(pos)
  
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
  
  cat("Generating sample covariance/correlation matrices...\n")
  result.Corr <- makeCorr(X, 1:N, h)
  Corr <- result.Corr$Corr
  if(fit.corr == TRUE) {
    sd.X <- result.Corr$sd.X
    cat("Use sample correlation matrices in model fitting\n")
  } else {
    sd.X <- rep(1, p)
    cat("Use sample covariance matrices in model fitting\n")
  }
  rm(result.Corr)
  
  cat("Estimating graphs...\n")
  
  Omega.list <- vector("list", K)
  Omega.rf.list <- vector("list", K)
  edge.num.list <- rep(0, K)
  edge.list <- vector("list", K)
  
  if(length(d) == 1 && d >= 0.5) {
    
    if(length(lambda) == 1) {
      
      result <- LGGM.global(pos, Corr, sd.X, fit.type, refit.type, lambda, epi.abs, epi.rel)
      
      cat("Complete all!\n")
      
    } else {
      
      lambda.list <- unique(lambda)
      
      for(lambda.l in lambda.list) {
        
        result.l <- LGGM.global(pos, Corr, sd.X, fit.type, refit.type, lambda.l, epi.abs, epi.rel)
        idx <- which(lambda == lambda.l)
        
        for(i in idx) {
          Omega.list[[i]] <- result.l$Omega.list[[i]]
          Omega.rf.list[[i]] <- result.l$Omega.rf.list[[i]]
          edge.num.list[i] <- result.l$edge.num.list[i]
          edge.list[[i]] <- result.l$edge.list[[i]]
        }
        
        cat("Complete: t =", round((pos[idx]-1) / (N-1), 2), "\n")
      }
      
      result <- new.env()
      result$Omega.list <- Omega.list
      result$Omega.rf.list <- Omega.rf.list
      result$edge.num.list <- edge.num.list
      result$edge.list <- edge.list
      result <- as.list(result)
    }
    
  } else{
    
    if(length(d) == 1) {
      d <- rep(d, K)
    }
    
    if(length(lambda) == 1) {
      lambda <- rep(lambda, K)
    }
    
    registerDoParallel(num.thread)
    
    result <- foreach(k = 1:K, .combine = "list", .multicombine = TRUE, .maxcombine = K, .export = c("LGGM.local")) %dopar%
      LGGM.local(pos[k], Corr, sd.X, fit.type, refit.type, d[k], lambda[k], epi.abs, epi.rel)
    
    for(k in 1:K) {
      
      result.k <- result[[k]]
      
      Omega.list[[k]] <- result.k$Omega
      Omega.rf.list[[k]] <- result.k$Omega.rf
      edge.num.list[[k]] <- result.k$edge.num
      edge.list[[k]] <- result.k$edge
    }
    
    result <- new.env()
    result$Omega.list <- Omega.list
    result$Omega.rf.list <- Omega.rf.list
    result$edge.num.list <- edge.num.list
    result$edge.list <- edge.list
    result <- as.list(result)
  }
  
  return(result)  
}