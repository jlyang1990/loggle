
# Graph estimation function for loggle #################################################################################
########################################################################################################################

# Input ###
# X: a p by N matrix containing list of observations
# pos: position of time points where graphs are estimated
# h: bandwidth in kernel function used to generate correlation matrices
# d: list of widths of neighborhood
# lambda: list of tuning parameters of Lasso penalty
# fit.type: 0: graphical Lasso estimation, 
#           1: pseudo likelihood estimation, 
#           2: sparse partial correlation estimation
# refit: whether to conduct model refitting
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion
# max.step: maximum steps in ADMM iteration
# detrend: whether to detrend variables in data matrix by subtracting kernel weighted moving average or overall average
# fit.corr: whether to use sample correlation matrix rather than sample covariance matrix in model fitting
# num.thread: number of threads
# print.detail: whether to print details in model fitting procedure

# Output ###
# Omega.list: if refit = TRUE: list of refitted precision matrices of length K (number of time points)
#             if refit = FALSE: list of estimated precision matrices of length K (number of time points)
# edge.num.list: list of detected edge numbers of length K
# edge.list: list of detected edges of length K

loggle <- function(X, pos = 1:ncol(X), h = 0.8*ncol(X)^(-1/5), d = 0.2, lambda = 0.25, fit.type = "pseudo", 
                   refit = TRUE, epi.abs = 1e-5, epi.rel = 1e-3, max.step = 500, detrend = TRUE, 
                   fit.corr = TRUE, num.thread = 1, print.detail = TRUE) {
  
  p <- dim(X)[1]
  N <- dim(X)[2]
  K <- length(pos)
  
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
  if(!length(d) %in% c(1, K)) {
    stop("d must be a scalar or a vector of the same length as 'pos'!")
  }
  if(!length(lambda) %in% c(1, K)) {
    stop("lambda must be a scalar or a vector of the same length as 'pos'!")
  }
  
  cat("Detrending each variable in data matrix...\n")
  X <- dataDetrend(X, detrend)
  
  cat("Generating sample covariance/correlation matrices...\n")
  result.Corr <- makeCorr(X, 1:N, h, fit.corr)
  Corr <- result.Corr$Corr
  sd.X <- result.Corr$sd.X
  rm(result.Corr)
  
  cat("Estimating graphs...\n")
  
  Omega.list <- vector("list", K)
  edge.num.list <- rep(0, K)
  edge.list <- vector("list", K)
  
  if(length(d) == 1 && d >= 0.5) {
    
    if(length(lambda) == 1) {
      
      result <- loggle.global(pos, Corr, sd.X, lambda, fit.type, refit, epi.abs, epi.rel, max.step)
      
      if(print.detail) {
        cat("Complete all!\n")
      }
      
    } else {
      
      lambda.list <- unique(lambda)
      
      for(lambda.l in lambda.list) {
        
        result.l <- loggle.global(pos, Corr, sd.X, lambda.l, fit.type, refit, epi.abs, epi.rel, max.step)
        idx <- which(lambda == lambda.l)
        
        for(i in idx) {
          Omega.list[[i]] <- result.l$Omega.list[[i]]
          edge.num.list[i] <- result.l$edge.num.list[i]
          edge.list[[i]] <- result.l$edge.list[[i]]
        }
        
        if(print.detail) {
          cat("Complete: t =", round((pos[idx]-1) / (N-1), 2), "\n")
        }
      }
      
      result <- list(Omega.list = Omega.list, edge.num.list = edge.num.list, edge.list = edge.list)
    }
    
  } else{
    
    if(length(d) == 1) {
      d <- rep(d, K)
    }
    
    if(length(lambda) == 1) {
      lambda <- rep(lambda, K)
    }
    
    if(num.thread > 1) {
      
      registerDoParallel(num.thread)
      
      result <- foreach(k=1:K, .combine="list", .multicombine=TRUE, .maxcombine=K, .export=c("loggle.local")) %dopar%
        loggle.local(pos[k], Corr, sd.X, d[k], lambda[k], fit.type, refit, epi.abs, epi.rel, max.step, print.detail)
      
    } else {
      
      result <- vector("list", K)
      for(k in 1:K) {
        result[[k]] <- loggle.local(pos[k], Corr, sd.X, d[k], lambda[k], fit.type, refit, epi.abs, epi.rel, max.step, 
                                    print.detail)
      }
    }
    
    for(k in 1:K) {
      
      result.k <- result[[k]]
      
      Omega.list[[k]] <- result.k$Omega
      edge.num.list[[k]] <- result.k$edge.num
      edge.list[[k]] <- result.k$edge
    }
    
    result <- list(Omega.list = Omega.list, edge.num.list = edge.num.list, edge.list = edge.list)
  }
  
  return(result)  
}


# Graph estimation function for local loggle ###########################################################################
########################################################################################################################

# Input ###
# pos: position of time point where graph is estimated
# Corr: list of kernel estimators of correlation matrices
# sd.X: list of standard deviations of variables
# d: width of neighborhood
# lambda: tuning parameter of Lasso penalty
# fit.type: 0: graphical Lasso estimation, 
#           1: pseudo likelihood estimation, 
#           2: sparse partial correlation estimation
# refit: whether to conduct model refitting
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion
# max.step: maximum steps in ADMM iteration
# print.detail: whether to print details in model fitting procedure

# Output ###
# Omega: if refit = TRUE: refitted precision matrix
#        if refit = FALSE: estimated precision matrix
# edge.num: detected edge number
# edge: detected edges

loggle.local <- function(pos, Corr, sd.X, d, lambda, fit.type, refit, epi.abs, epi.rel, max.step, print.detail) {
  
  p <- dim(Corr)[1]
  N <- dim(Corr)[3]
  
  Nd.index <- max(1, ceiling(((pos-1)/(N-1)-d)*(N-1)-1e-5)+1) : min(N, floor(((pos-1)/(N-1)+d)*(N-1)+1e-5)+1)
  Nd <- length(Nd.index)
  Nd.pos <- which(Nd.index == pos)
  Nd.pos.c <- Nd.pos - 1
  Nd.pos.l <- 1
  
  Corr.sq <- rowSums(Corr[, , Nd.index, drop = FALSE]^2, dims = 2)
  
  Z.vec <- rep(0, p*p)
  
  lambda <- sqrt(Nd) * lambda
  rho <- lambda
  
  #detect block diagonal structure
  adj.mat <- (Corr.sq > lambda^2)
  diag(adj.mat) <- 1
  graph <- graph.adjacency(adj.mat)
  cluster <- clusters(graph)
  member <- cluster$membership
  csize <- cluster$csize
  no <- cluster$no
  member.index <- sort(member, index.return = T)$ix - 1
  csize.index <- c(0, cumsum(csize))
  
  result <- .C("ADMM_simple",
               as.double(Corr[, , Nd.index]),
               Z.vec = as.double(Z.vec),
               as.integer(p),
               as.integer(Nd),
               as.integer(Nd.pos.c),
               as.integer(Nd.pos.l),
               as.integer(member.index),
               as.integer(csize.index),
               as.integer(no),
               as.double(lambda),
               as.integer(fit.type),
               as.double(rho),
               as.double(epi.abs),
               as.double(epi.rel),
               as.integer(max.step)
  )
  
  Omega <- Matrix(result$Z.vec, p, p, sparse = T)
  
  edge <- which(as.matrix(Omega) != 0, arr.ind = T)
  edge <- edge[(edge[, 1] - edge[, 2]) > 0, , drop = F]
  edge.num <- nrow(edge)
  
  if(refit) {
    
    edge.zero <- which(as.matrix(Omega) == 0, arr.ind = T)
    edge.zero <- edge.zero[(edge.zero[, 1] - edge.zero[, 2]) > 0, , drop = F]
    if(nrow(edge.zero) == 0) {
      edge.zero = NULL
    }
    
    Sigma <- diag(sd.X) %*% Corr[, , pos] %*% diag(sd.X)
    Omega <- glasso::glasso(s = Sigma, rho = 1e-10, zero = edge.zero)$wi  # set rho = 1e-10 instead of 0 to avoid warning
    if(det(Omega) < 0) {
      Omega <- glasso::glasso(s = Sigma, rho = 1e-10, zero = edge.zero, thr = 5*1e-5)$wi
    }
    Omega <- Matrix(Omega, sparse = T)
  }
  
  if(print.detail) {
    cat("Complete: t =", round((pos-1) / (N-1), 2), "\n")
  }
  
  result <- list(Omega = Omega, edge.num = edge.num, edge = edge)
  return(result)
}


# Graph estimation function for global loggle ##########################################################################
########################################################################################################################

# Input ###
# pos: position of time points where graphs are estimated
# Corr: list of kernel estimators of correlation matrices
# sd.X: list of standard deviations of variables
# d: width of neighborhood
# lambda: tuning parameter of Lasso penalty
# fit.type: 0: graphical Lasso estimation, 
#           1: pseudo likelihood estimation, 
#           2: sparse partial correlation estimation
# refit: whether to conduct model refitting
# epi.abs: absolute tolerance in ADMM stopping criterion
# epi.rel: relative tolerance in ADMM stopping criterion
# max.step: maximum steps in ADMM iteration

# Output ###
# Omega: if refit = TRUE: list of refitted precision matrices of length K (number of time points)
#        if refit = FALSE: list of estimated precision matrices of length K
# edge.num: list of detected edge numbers of length K
# edge: list of detected edges of length K

loggle.global <- function(pos, Corr, sd.X, lambda, fit.type, refit, epi.abs, epi.rel, max.step) {
  
  p <- dim(Corr)[1]
  N <- dim(Corr)[3]
  K <- length(pos)
  
  N.index.c <- 0:(N-1)
  pos.c <- pos - 1
  
  Corr.sq <- rowSums(Corr^2, dims = 2)
  
  Z.vec <- rep(0, p*p*K)
  
  lambda <- sqrt(N) * lambda
  rho <- lambda
  
  #detect block diagonal structure
  adj.mat <- (Corr.sq > lambda^2)
  diag(adj.mat) <- 1
  graph <- graph.adjacency(adj.mat)
  cluster <- clusters(graph)
  member <- cluster$membership
  csize <- cluster$csize
  no <- cluster$no
  member.index <- sort(member, index.return = T)$ix - 1
  csize.index <- c(0, cumsum(csize))
  
  result <- .C("ADMM_simple",
               as.double(Corr),
               Z.vec = as.double(Z.vec),
               as.integer(p),
               as.integer(N),
               as.integer(pos.c),
               as.integer(K),
               as.integer(member.index),
               as.integer(csize.index),
               as.integer(no),
               as.double(lambda),
               as.integer(fit.type),
               as.double(rho),
               as.double(epi.abs),
               as.double(epi.rel),
               as.integer(max.step)
  )
  
  Omega.list <- sapply(1:K, function(k) Matrix(result$Z.vec[(p*p*(k-1) + 1) : (p*p*k)], p, p, sparse = T))
  
  edge <- which(as.matrix(Omega.list[[1]]) != 0, arr.ind = T)
  edge <- edge[(edge[, 1] - edge[, 2]) > 0, , drop = F]
  edge.list <- rep(list(edge), K)
  edge.num.list <- rep(nrow(edge), K)
  
  if(refit) {
    
    edge.zero <- which(as.matrix(Omega.list[[1]]) == 0, arr.ind = T)
    edge.zero <- edge.zero[(edge.zero[, 1] - edge.zero[, 2]) > 0, , drop = F]
    if(nrow(edge.zero) == 0) {
      edge.zero = NULL
    }
    
    Omega.list <- vector("list", K)
    for(k in 1:K) {
      Sigma <- diag(sd.X) %*% Corr[, , pos[k]] %*% diag(sd.X)
      Omega.list[[k]] <- glasso::glasso(s = Sigma, rho = 1e-10, zero = edge.zero)$wi
      if(det(Omega.list[[k]]) < 0) {
        Omega.list[[k]] <- glasso::glasso(s = Sigma, rho = 1e-10, zero = edge.zero, thr = 5*1e-5)$wi
      }
      Omega.list[[k]] <- Matrix(Omega.list[[k]], sparse = T)
    }
  }
  
  result <- list(Omega.list = Omega.list, edge.num.list = edge.num.list, edge.list = edge.list)
  return(result)
}


# Data detrending function #############################################################################################
########################################################################################################################

# Input ###
# X: a p by N matrix containing list of observations
# detrend: whether to detrend variables in data matrix by subtracting kernel weighted moving average or overall average 

# Output ###
# X: a p by N matrix containing list of observations after detrending

dataDetrend <- function(X, detrend) {
  
  p <- nrow(X)
  N <- ncol(X)
  
  for(i in 1:p) {
    if(detrend) {
      X[i, ] <- X[i, ] - sm::sm.regression(1:N, X[i, ], ngrid = N, display = "none")$estimate
    } else {
      X[i, ] <- X[i, ] - mean(X[i, ])
    }
  }
  
  return(X)
}


# Covariance/Correlation matrix generation function ####################################################################
########################################################################################################################

# Input ###
# X: a p by N matrix containing list of observations
# pos: position of observations used to generate correlation matrices
# h: bandwidth in kernel function used to generate correlation matrices
# fit.corr: whether to use sample correlation matrix rather than sample covariance matrix in model fitting

# Output ###
# Corr: list of kernel estimators of correlation matrices
# sd.X: list of standard deviations of variables when fit.corr = TRUE, list of 1's when fit.corr = FALSE

makeCorr <- function(X, pos, h, fit.corr) {
  
  p <- nrow(X)
  N <- ncol(X)
  L <- length(pos)
  
  sd.X <- rep(1, p)
  
  if(fit.corr) {
    for(i in 1:p) {
      sd.X[i] <- sd(X[i, ])
      X[i, ] <- X[i, ] / sd(X[i, ])
    }
  }
  
  Corr <- array(0, c(p, p, N))
  
  for(i in 1:N) {
    Kh <- pmax(3/4 * (1 - ((pos - i) / ((N - 1) * h))^2), 0)
    omega <- Kh / sum(Kh)
    index <- which(omega != 0)
    X_pos <- X[, pos[index]]
    Corr[, , i] <- (X_pos * rep(omega[index], rep(p, length(index)))) %*% t(X_pos)
  }
  
  result <- list(Corr = Corr, sd.X = sd.X)
  return(result)
}