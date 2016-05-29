
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

LGGM.local = function(pos, Corr, Sigma, d, lambda, pseudo.fit, pseudo.refit, epi.abs, epi.rel){
  
  p = dim(Sigma)[1]; N = dim(Sigma)[3]
    
  Nd.index = max(1, ceiling(((pos-1)/(N-1)-d)*(N-1)-1e-5)+1):min(N, floor(((pos-1)/(N-1)+d)*(N-1)+1e-5)+1)
  Nd.index.c = Nd.index - 1; Nd = length(Nd.index)
  Nd.pos = which(Nd.index == pos); Nd.pos.c = Nd.pos - 1; Nd.pos.l = 1
    
  Corr.sq = apply(Corr^2, c(1, 2), sum)
    
  Z.ini.vec = rep(0, p*p); Z.pos.vec = rep(0, p*p); U.vec = rep(0, p*p); edge = 0
    
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
              as.double(Corr),
              as.double(Sigma),
              Z.ini.vec = as.double(Z.ini.vec),
              Z.pos.vec = as.double(Z.pos.vec),
              as.double(U.vec),
              as.double(lambda),
              as.double(rho),
              as.double(epi.abs),
              as.double(epi.rel),
              as.integer(pseudo.fit),
              as.integer(pseudo.refit),
              as.integer(edge)
  )
      
  Omega.lg = Matrix(result$Z.ini.vec, p, p, sparse=T)
  Omega.rf = Matrix(result$Z.pos.vec, p, p, sparse=T)
  edge = result$edge
  S = which(Omega.rf!=0, arr.ind=T); S = S[(S[, 1] - S[, 2])>0, , drop=F]
  cv = sum(c(t(Sigma[, , pos]))*c(Omega.rf)) - log(det(Omega.rf))
    
  cat("complete: t =", round((pos-1)/(N-1), 2), "\n")
  
  return(list(Omega.lg, Omega.rf, edge, S, cv))
}


#simulation study for global group-lasso (d=1)################################################################################################################################################################################################################

LGGM.global = function(pos, Corr, Sigma, lambda, pseudo.fit, pseudo.refit, epi.abs, epi.rel){
  
  p = dim(Sigma)[1]; N = dim(Sigma)[3]
  
  N.index.c = 0:(N-1)
  pos.c = pos-1
  
  Corr.sq = apply(Corr[, , N.index]^2, c(1, 2), sum)
  
  Z.ini.vec = rep(0, p*p*K); Z.pos.vec = rep(0, p*p*K); U.vec = rep(0, p*p*K); edge = 0
  
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
              Z.ini.vec = as.double(Z.ini.vec),
              Z.pos.vec = as.double(Z.pos.vec),
              as.double(U.vec),
              as.double(lambda),
              as.double(rho),
              as.double(epi.abs),
              as.double(epi.rel),
              as.integer(pseudo.fit),
              as.integer(pseudo.refit),
              as.integer(edge)
  )
  
  Omega.lg.list = sapply(1:K, function(k) Matrix(result$Z.ini.vec[(p*p*(k-1)+1):(p*p*k)], p, p, sparse=T))
  Omega.rf.list = sapply(1:K, function(k) Matrix(result$Z.pos.vec[(p*p*(k-1)+1):(p*p*k)], p, p, sparse=T))
  edge.list = rep(result$edge, K)
  S = which(Omega.rf!=0, arr.ind=T); S = S[(S[, 1] - S[, 2])>0, , drop=F]; S.list = rep(list(S), K)
  cv.list = sapply(1:K, function(k) sum(c(t(Sigma[, , pos[k]]))*c(Omega.rf)) - log(det(Omega.rf)))
  
  return(list(Omega.lg.list, Omega.rf.list, edge.list, S.list, cv.list))
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

LGGM = function(X, pos, pseudo.fit, h, d, lambda, epi.abs, epi.rel, corr.ind = TRUE, pseudo.refit, thres){
  
  p = dim(X)[1]; N = dim(X)[2]; K = length(pos)
  
  Sigma = gene.Sigma(X, pos, h)
  if(corr.ind == TRUE){
    Corr = gene.corr(X, pos, h)
  }else{
    Corr = Sigma
  }
  
  Nd.index = max(1, ceiling(((pos-1)/(N-1)-d)*(N-1)-1e-5)+1):min(N, floor(((pos-1)/(N-1)+d)*(N-1)+1e-5)+1)
  Corr = Corr[, , Nd.index]
  Sigma = Sigma[, , Nd.index]
  
  if(d!=1){
    
    result = foreach(k = 1:K, .combine = "list", .multicombine = TRUE, .maxcombine = K, .export = c("LGGM.pos")) %dopar%
      LGGM.local(pos[k], Corr, Sigma, d, lambda, pseudo.fit, pseudo.refit, epi.abs, epi.rel)
    
    S.list = vector("list", K)
    Omega.lg.list = vector("list", K); Omega.rf.list = vector("list", K)
    edge.list = rep(0, K); cv.list = rep(0, K)
    
    for(k in 1:K){
      
      result.k = result[[k]]
      
      Omega.lg.list[k] = result.k[[1]]
      Omega.rf.list[k] = result.k[[2]]
      edge.list[k] = result.k[[3]]
      S.list[k] = result.k[[4]]
      cv.list[k] = result.k[[5]]
    }
  
  }else{
    
    result = LGGM.global(pos, Corr, Sigma, lambda, pseudo.fit, pseudo.refit, epi.abs, epi.rel)
    
    Omega.lg.list = result[[1]]
    Omega.rf.list = result[[2]]
    edge.list = result[[3]]
    S.list = result[[4]]
    cv.list = result[[5]]
  }
  
  result = new.env()
  result$Omega.lg.list = Omega.lg.list
  result$Omega.rf.list = Omega.rf.list
  result$edge.list = edge.list
  result$S.list = S.list
  result$cv.list = cv.list
  
  result = as.list(result)
  
  return(result)  
}