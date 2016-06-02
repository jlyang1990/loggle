
#Sigma generation function####################################################################################################################################################################################################################################

#Input###
#X: list of observations
#pos: position of observations used to generate Sigma
#h: bandwidth in kernel function

#Output###
#Sigma: list of kernel estimators of covariance matrices

gene.Sigma = function(X, pos, h){
  
  p = nrow(X); N = ncol(X); L = length(pos)
  
  Sigma = array(0, c(p, p, N))
  
  for(i in 1:N){
    Kh = pmax(3/4*(1-((pos-i)/((N-1)*h))^2), 0); omega = Kh/sum(Kh); index = which(omega!=0)
    Sigma[, , i] = X[, pos[index]]%*%diag(omega[index])%*%t(X[, pos[index]])
  }
  
  return(Sigma)
}


#Correlation matrix generation function#######################################################################################################################################################################################################################

gene.corr= function(X, pos, h){
  
  p = nrow(X); N = ncol(X); L = length(pos)
  
  for(i in 1:p){X[i, ] = scale(X[i, ])}
  
  Corr = array(0, c(p, p, N))
  
  for(i in 1:N){
    Kh = pmax(3/4*(1-((pos-i)/((N-1)*h))^2), 0); omega = Kh/sum(Kh); index = which(omega!=0)
    Corr[, , i] = X[, pos[index]]%*%diag(omega[index])%*%t(X[, pos[index]])
  }
  
  return(Corr)
}


#simulation study for local group-lasso at time point pos#####################################################################################################################################################################################################

#Input###
#pos: position of time point where graph is estimated
#d.l: list of d's
#lambda.c: list of lambda's
#epi.abs, epi.rel: constants in ADMM stopping criterion
#pseudo.fit: 0: local group graphical lasso, 1: pseudo-likelihood group lasso (asymmetry), 2: pseudo-likelihood group lasso (symmetry), 3: SPACE
#pseudo.refit: 0: likelihood refit, 1: pseudo-likelihood refit
#thres: grid search stops when number of detected edges larger than thres*p

#Output###
#S.list: D (number of d's) by L (number of lambda's) lists of edges
#Omega.lg.list: D by L lists of estimated precision matrices
#Omega.rf.list: D by L lists of refitted precision matrices
#edge: D by L numbers of edges
#time: the total time cost
#record.list: D by L time costs in C
#time.list: D time costs w.r.t. D d's in R

LGGM.local = function(pos, Corr, Sigma, d, lambda, fit.type, refit.type, epi.abs, epi.rel){
  
  p = dim(Sigma)[1]; N = dim(Sigma)[3]
  
  Nd.index = max(1, ceiling(((pos-1)/(N-1)-d)*(N-1)-1e-5)+1):min(N, floor(((pos-1)/(N-1)+d)*(N-1)+1e-5)+1)
  Nd.index.c = Nd.index-1; Nd = length(Nd.index)
  Nd.pos = which(Nd.index == pos); Nd.pos.c = Nd.pos-1; Nd.pos.l = 1
    
  Corr.sq = apply(Corr[, , Nd.index]^2, c(1, 2), sum)
    
  Z.vec = rep(0, p*p*Nd); Z.pos.vec = rep(0, p*p); U.vec = rep(0, p*p*Nd); edge.num = 0
    
  lambda = sqrt(Nd)*lambda; rho = lambda
    
  #detect block diagonal structure
  adj.mat = (Corr.sq>lambda^2); diag(adj.mat) = 1; graph = graph.adjacency(adj.mat)
  cluster = clusters(graph); member = cluster$membership; csize = cluster$csize; no = cluster$no
  member.index = sort(member, index.return=T)$ix-1; csize.index = c(0, cumsum(csize))
    
  result = .C("ADMM_cluster",
              as.integer(p),
              as.integer(member.index),
              as.integer(csize.index),
              as.integer(no),
              as.integer(Nd),
              as.integer(Nd.pos.c),
              as.integer(Nd.pos.l),
              as.double(Corr[, , Nd.index]),
              as.double(Sigma[, , Nd.index]),
              Z.vec = as.double(Z.vec),
              Z.pos.vec = as.double(Z.pos.vec),
              as.double(U.vec),
              as.double(lambda),
              as.double(rho),
              as.double(epi.abs),
              as.double(epi.rel),
              as.integer(fit.type),
              as.integer(refit.type),
              edge.num = as.integer(edge.num)
  )
      
  Omega = Matrix(result$Z.vec[(p*p*(Nd.pos-1)+1):(p*p*Nd.pos)], p, p, sparse=T)
  Omega.rf = Matrix(result$Z.pos.vec, p, p, sparse=T)
  edge.num = result$edge.num
  edge = which(Omega.rf!=0, arr.ind=T); edge = edge[(edge[, 1] - edge[, 2])>0, , drop=F]
  cv = sum(c(t(Sigma[, , pos]))*c(matrix(Omega.rf))) - log(det(Omega.rf))
    
  cat("complete: t =", round((pos-1)/(N-1), 2), "\n")
  
  result = new.env()
  result$Omega = Omega
  result$Omega.rf = Omega.rf
  result$edge.num = edge.num
  result$edge = edge
  result$cv = cv
  
  result = as.list(result)
  
  return(result)
}


#simulation study for global group-lasso (d=1)################################################################################################################################################################################################################

LGGM.global = function(pos, Corr, Sigma, lambda, fit.type, refit.type, epi.abs, epi.rel){
  
  p = dim(Sigma)[1]; N = dim(Sigma)[3]; K = length(pos)
  
  N.index.c = 0:(N-1)
  pos.c = pos-1
  
  Corr.sq = apply(Corr^2, c(1, 2), sum)
  
  Z.vec = rep(0, p*p*N); Z.pos.vec = rep(0, p*p*K); U.vec = rep(0, p*p*N); edge.num = 0
  
  lambda = sqrt(N)*lambda; rho = lambda
  
  #detect block diagonal structure
  adj.mat = (Corr.sq>lambda^2); diag(adj.mat) = 1; graph = graph.adjacency(adj.mat)
  cluster = clusters(graph); member = cluster$membership; csize = cluster$csize; no = cluster$no
  member.index = sort(member, index.return=T)$ix-1; csize.index = c(0, cumsum(csize))
  
  result = .C("ADMM_cluster", 
              as.integer(p),
              as.integer(member.index),
              as.integer(csize.index),
              as.integer(no),
              as.integer(N),
              as.integer(pos.c),
              as.integer(K),
              as.double(Corr),
              as.double(Sigma),
              Z.vec = as.double(Z.vec),
              Z.pos.vec = as.double(Z.pos.vec),
              as.double(U.vec),
              as.double(lambda),
              as.double(rho),
              as.double(epi.abs),
              as.double(epi.rel),
              as.integer(fit.type),
              as.integer(refit.type),
              edge.num = as.integer(edge.num)
  )
  
  Z.vec = array(result$Z.vec, c(p, p, N))[, , pos]
  Omega.list = sapply(1:K, function(k) Matrix(Z.vec[, , k], p, p, sparse=T))
  Omega.rf.list = sapply(1:K, function(k) Matrix(result$Z.pos.vec[(p*p*(k-1)+1):(p*p*k)], p, p, sparse=T))
  edge.num.list = rep(result$edge.num, K)
  edge = which(Omega.rf.list[[1]]!=0, arr.ind=T); edge = edge[(edge[, 1] - edge[, 2])>0, , drop=F]; edge.list = rep(list(edge), K)
  cv = mean(sapply(1:K, function(k) sum(c(t(Sigma[, , pos[k]]))*c(matrix(Omega.rf.list[[k]]))) - log(det(Omega.rf.list[[k]]))))
  
  result = new.env()
  result$Omega.list = Omega.list
  result$Omega.rf.list = Omega.rf.list
  result$edge.num.list = edge.num.list
  result$edge.list = edge.list
  result$cv = cv
  
  result = as.list(result)
  
  return(result)
}


#simulation study for local group-lasso (main function)#######################################################################################################################################################################################################

#Input###
#pos: list of positions of time points where graphs are estimated
#h: bandwidth in kernel function
#X: list of observations
#corr.ind: 0: use Sigma in model fitting, 1: use Corr in model fitting
#d.l: list of d's
#lambda.c: list of lambda's
#epi.abs, epi.rel: constants in ADMM stopping criterion
#pseudo.fit: 0: local group graphical lasso, 1: pseudo-likelihood group lasso (asymmetry), 2: pseudo-likelihood group lasso (symmetry), 3: SPACE
#pseudo.refit: 0: likelihood refit, 1: pseudo-likelihood refit
#thres: grid search stops when number of detected edges larger than thres*p

#Output###
#S.list: D (number of d's) by L (number of lambda's) by K (number of time points) lists of edges
#Omega.lg.list: D by L by K lists of estimated precision matrices
#Omega.rf.list: D by L by K lists of refitted precision matrices
#edge: D by L by K numbers of edges
#time: K time costs w.r.t. K time points in R
#record.list: D by L by K time costs in C
#time.list: D by K time costs in R

LGGM = function(X, pos = 1:ncol(X), fit.type = "glasso", refit.type = "glasso", h = 0.2, d, lambda, epi.abs = 1e-5, epi.rel = 1e-3, corr = TRUE, num.core = 1){
  
  p = dim(X)[1]; N = dim(X)[2]; K = length(pos)
  
  if(fit.type == "glasso"){
    fit.type = 0
  }
  if(fit.type == "pseudo"){
    fit.type = 1
  }
  if(fit.type == "space"){
    fit.type = 2
  }
  
  if(refit.type == "glasso"){
    refit.type = 0
  }
  if(refit.type == "pseudo"){
    refit.type = 1
  }
  
  Sigma = gene.Sigma(X, 1:N, h)
  if(corr == TRUE){
    Corr = gene.corr(X, 1:N, h)
  }else{
    Corr = Sigma
  }
  
  if(length(d) == 1 && d>=0.5){
    
    result = LGGM.global(pos, Corr, Sigma, lambda, fit.type, refit.type, epi.abs, epi.rel)
    
  } else{
    
    if(length(d) == 1){
      d = rep(d, K)
    }
    
    if(length(lambda) == 1){
      lambda = rep(lambda, K)
    }
    
    registerDoParallel(num.core)
    
    result = foreach(k = 1:K, .combine = "list", .multicombine = TRUE, .maxcombine = K, .export = c("LGGM.local")) %dopar%
      LGGM.local(pos[k], Corr, Sigma, d[k], lambda[k], fit.type, refit.type, epi.abs, epi.rel)
    
    stopImplicitCluster()
    
    Omega.list = vector("list", K); Omega.rf.list = vector("list", K)
    edge.num.list = vector("list", K); edge.list = rep(0, K); cv.list = rep(0, K)
    
    for(k in 1:K){
      
      result.k = result[[k]]
      
      Omega.list[[k]] = result.k$Omega
      Omega.rf.list[[k]] = result.k$Omega.rf
      edge.num.list[[k]] = result.k$edge.num
      edge.list[[k]] = result.k$edge
      cv.list[[k]] = result.k$cv
    }
    
    cv = mean(cv.list)
    
    result = new.env()
    result$Omega.list = Omega.list
    result$Omega.rf.list = Omega.rf.list
    result$edge.num.list = edge.num.list
    result$edge.list = edge.list
    result$cv = cv
    
    result = as.list(result)
  }
  
  return(result)  
}