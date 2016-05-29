
library(MASS); library(Matrix); library(sparseMVN); library(matrixcalc); library(foreach); library(doParallel); library(igraph)

dyn.load("ADMM.so")

#list function################################################################################################################################################################################################################################################

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}


#Data generation function#####################################################################################################################################################################################################################################

#Input###
#p: dimension of each observation
#N: number of observations (number of time points in observations)
#alpha: value of threshold in soft-thresholding

#Output###
#Omega: list of true precision matrices
#X: list of observations
#edge: list of numbers of edges in true precision matrices

gene.data = function(p, N, alpha){
  
  pd = F
  
  while(pd==F){
    
    t = seq(0, 1, length=N)
    B1 = matrix(rnorm(p^2, 0, 1/2), p, p); B1[upper.tri(B1)] = 0
    B2 = matrix(rnorm(p^2, 0, 1/2), p, p); B2[upper.tri(B2)] = 0
    B3 = matrix(rnorm(p^2, 0, 1/2), p, p); B3[upper.tri(B3)] = 0
    B4 = matrix(rnorm(p^2, 0, 1/2), p, p); B4[upper.tri(B4)] = 0
    Omega = array(0, c(p, p, N)); edge = rep(0, N); X = matrix(0, p, N)
    uni = runif(4, -0.5, 0.5)
    
    for(i in 1:N){
      
      G = (B1*sin(pi*t[i]/2+uni[1]) + B2*cos(pi*t[i]/2+uni[2]) + B3*sin(pi*t[i]/4+uni[3]) + B4*cos(pi*t[i]/4+uni[4]))/2
      Omega.o = G%*%t(G)
      Omega.o = Omega.o%*%diag(1/sqrt((1:p))); Omega.o[upper.tri(Omega.o)] = 0
      Omega.o = Omega.o + t(Omega.o); diag(Omega.o) = diag(Omega.o)/(2*sqrt(1:p)) + log(p, 10)/4
      Omega.temp1 = matrix(1, p, p) - alpha/abs(Omega.o)
      Omega.temp2 = matrix(1, p, p) - alpha/(2*abs(Omega.o))
      Omega.n = Omega.temp2*(Omega.temp1>0)*Omega.o
      diag(Omega.n) = diag(Omega.o)
      
      if(is.positive.definite(Omega.n)==F){pd = F; break}
      
      Omega[, , i] = Omega.n
      Omega.n = Matrix(Omega.n, sparse=T)
      edge[i] = (sum(!(Omega.n==0))-p)/2
      X[, i] = rmvn.sparse(1, rep(0, p), Cholesky(Omega.n))# + 0*rnorm(p, sd = 0.2)
      if(i==N){pd = T}
    }
  }
  return(list(Omega, X, edge))
}


#Data generation function (time-invariant)####################################################################################################################################################################################################################

gene.data.2 = function(p, N){
  pd = F
  while(pd==F){
    t = seq(0, 1, length=N)
    B = matrix((runif(p^2)<(2/p)), p, p); B[upper.tri(B)] = 0; B = B + t(B); diag(B) = 1
    C = matrix(runif(p^2), p, p); C[upper.tri(C)] = 0; C = C + t(C); diag(C) = diag(C)/2
    Omega = array(0, c(p, p, N)); edge = rep(0, N); X = matrix(0, p, N)
    for(i in 1:N){
      Omega.n = sin(2*pi*t[i] + C)*B
      diag(Omega.n) = abs(diag(Omega.n)) +log(p, 10)
      if(is.positive.definite(Omega.n)==F){pd = F; break}
      Omega[, , i] = Omega.n
      Omega.n = Matrix(Omega.n, sparse=T)
      edge[i] = (sum(!(Omega.n==0))-p)/2
      X[, i] = rmvn.sparse(1, rep(0, p), Cholesky(Omega.n))
      if(i==N){pd = T}
    }
  }
  return(list(Omega, X, edge))
}


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

simu.lglasso.pos = function(pos, Corr, Sigma, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres){
  
  ptm = proc.time()
  
  p = dim(Sigma)[1]; N = dim(Sigma)[3]; D = length(d.l); L = length(lambda.c)
  
  S.list = matrix(vector("list", 1), L, D)
  Omega.lg.list = matrix(vector("list", 1), L, D); Omega.rf.list = matrix(vector("list", 1), L, D)
  edge = matrix(0, L, D)
  
  time.list = rep(0, D); record.list = matrix(0, L, D)
  
  for(j in 1:D){
    
    d = d.l[j]; epi.abs.d = epi.abs[j]; epi.rel.d = epi.rel[j]
    
    Nd.index = max(1, ceiling(((pos-1)/(N-1)-d)*(N-1)-1e-5)+1):min(N, floor(((pos-1)/(N-1)+d)*(N-1)+1e-5)+1)
    Nd.index.c = Nd.index - 1; Nd = length(Nd.index); Nd.pos = which(Nd.index == pos)
    Nd.pos.c = Nd.pos - 1; Nd.pos.l = 1
    
    Corr.sq = apply(Corr[, , Nd.index]^2, c(1, 2), sum)
    
    Z.ini.vec = rep(0, p*p*L); Z.pos.vec = rep(0, p*p*L)
    
    lambda = sqrt(Nd)*lambda.c; rho = lambda
    
    member.index.list = rep(0, p*L); no.list = rep(0, L); csize.index.list = c()
    
    diff2 = rep(0, L)
    
    #detect block diagonal structure
    for(l in L:1){
      
      adj.mat = (Corr.sq>lambda[l]^2); diag(adj.mat) = 1; graph = graph.adjacency(adj.mat)
      cluster = clusters(graph); member = cluster$membership; csize = cluster$csize; no = cluster$no
      member.index = sort(member, index.return=T)$ix-1; csize.index = c(0, cumsum(csize))
      
      member.index.list[(p*(l-1)+1):(p*l)] = member.index
      no.list[l] = no
      csize.index.list = c(csize.index.list, csize.index)
    }
      
    ptm1 = proc.time()
      
    result = .C("ADMM_lambda", 
                  as.integer(p),
                  as.integer(member.index.list),
                  as.integer(csize.index.list),
                  as.integer(no.list),
                  as.integer(Nd),
                  as.integer(Nd.pos.c),
                  as.integer(Nd.pos.l),
                  as.double(Corr[, , Nd.index]),
                  as.double(Sigma[, , Nd.index]),
                  Z.ini.vec = as.double(Z.ini.vec),
                  Z.pos.vec = as.double(Z.pos.vec),
                  as.integer(L),
                  as.double(lambda),
                  as.double(rho),
                  as.double(epi.abs.d),
                  as.double(epi.rel.d),
                  as.integer(pseudo.fit),
                  as.integer(pseudo.refit),
                  as.double(thres),
                  diff2 = as.double(diff2)
    )
      
    diff1 = (proc.time() - ptm1)[3]
      
    Z.ini.vec = result$Z.ini.vec
    Z.pos.vec = result$Z.pos.vec
      
    for(l in L:1){
      
      Omega.lg = matrix(Z.ini.vec[(p*p*(l-1)+1):(p*p*l)], p, p)
      Omega.rf = matrix(Z.pos.vec[(p*p*(l-1)+1):(p*p*l)], p, p)
        
      S = which(Omega.rf!=0, arr.ind=T); S = S[(S[, 1] - S[, 2])>0, , drop=F]; S.L = nrow(S)
        
      S.list[[l, j]] = S
      Omega.lg.list[[l, j]] = Matrix(Omega.lg, sparse = T)
      Omega.rf.list[[l, j]] = Matrix(Omega.rf, sparse = T)
      edge[l, j] = S.L
    }
      
    time.list[j] = diff1; record.list[, j] = result$diff2
    
    cat("complete: d =", d, "t =", round((pos-1)/(N-1), 2), "\n")
  }
  
  time = (proc.time() - ptm)[3]
  
  return(list(S.list, Omega.lg.list, Omega.rf.list, edge, time, record.list, time.list))
}


#simulation study for global group-lasso (d=1)################################################################################################################################################################################################################

simu.gglasso = function(pos, Corr, Sigma, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres){
  
  p = dim(Sigma)[1]; N = dim(Sigma)[3]; K = length(pos); L = length(lambda.c)
  
  S.list = matrix(vector("list", 1), L, K)
  Omega.lg.list = matrix(vector("list", 1), L, K); Omega.rf.list = matrix(vector("list", 1), L, K)
  edge = matrix(0, L, K)
  
  N.index.c = 0:(N-1)
  pos.c = pos-1
  
  Corr.sq = apply(Corr^2, c(1, 2), sum)
  
  Z.ini.vec = rep(0, p*p*K*L); Z.pos.vec = rep(0, p*p*K*L)
  
  lambda = sqrt(N)*lambda.c; rho = lambda
  
  member.index.list = rep(0, p*L); no.list = rep(0, L); csize.index.list = c()
  
  diff2 = rep(0, L)
  
  for(l in L:1){
    
    adj.mat = (Corr.sq>lambda[l]^2); diag(adj.mat) = 1; graph = graph.adjacency(adj.mat)
    cluster = clusters(graph); member = cluster$membership; csize = cluster$csize; no = cluster$no
    member.index = sort(member, index.return=T)$ix-1; csize.index = c(0, cumsum(csize))
    
    member.index.list[(p*(l-1)+1):(p*l)] = member.index
    no.list[l] = no
    csize.index.list = c(csize.index.list, csize.index)
  }
  
  result = .C("ADMM_lambda", 
              as.integer(p),
              as.integer(member.index.list),
              as.integer(csize.index.list),
              as.integer(no.list),
              as.integer(N),
              as.integer(pos.c),
              as.integer(K),
              as.double(Corr),
              as.double(Sigma),
              Z.ini.vec = as.double(Z.ini.vec),
              Z.pos.vec = as.double(Z.pos.vec),
              as.integer(L),
              as.double(lambda),
              as.double(rho),
              as.double(epi.abs),
              as.double(epi.rel),
              as.integer(pseudo.fit),
              as.integer(pseudo.refit),
              as.double(thres),
              diff2 = as.double(diff2)
  )
  
  Z.ini.vec = result$Z.ini.vec
  Z.pos.vec = result$Z.pos.vec
  
  for(l in L:1){
    
    Omega.gg = array(Z.ini.vec[(p*p*K*(l-1)+1):(p*p*K*l)], c(p, p, K))
    Omega.rf = array(Z.pos.vec[(p*p*K*(l-1)+1):(p*p*K*l)], c(p, p, K))
    
    S = which(Omega.rf[, , 1]!=0, arr.ind=T); S = S[(S[, 1] - S[, 2])>0, , drop=F]; S.L = nrow(S)
  
    edge[l, ] = S.L
    
    for(k in 1:K){
      
      S.list[[l, k]] = S
      Omega.lg.list[[l, k]] = Matrix(Omega.gg[, , k], sparse = T)
      Omega.rf.list[[l, k]] = Matrix(Omega.rf[, , k], sparse = T)
    }
  }
  
  cat("complete: d = 1", "\n")
  
  return(list(S.list, Omega.lg.list, Omega.rf.list, edge))
}


#simulation study for local group-lasso (main function)#######################################################################################################################################################################################################

#Input###
#pos: list of positions of time points where graphs are estimated
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

simu.lglasso = function(pos, Corr, Sigma, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres){
  
  p = dim(Sigma)[1]; N = dim(Sigma)[3]; K = length(pos); D = length(d.l); L = length(lambda.c)
  
  S.list = array(vector("list", 1), c(L, D, K))
  Omega.lg.list = array(vector("list", 1), c(L, D, K)); Omega.rf.list = array(vector("list", 1), c(L, D, K))
  edge = array(0, c(L, D, K)); time = rep(0, K)
  
  if(d.l[D]==1){
    
    record.list = array(0, c(L, D-1, K)); time.list = matrix(0, D-1, K)
    
    result.lg = foreach(k = 1:K, .combine = "list", .multicombine = TRUE, .maxcombine = K, .export = c("simu.lglasso.pos")) %dopar%
      simu.lglasso.pos(pos[k], Corr, Sigma, d.l[-D], lambda.c, epi.abs[-D], epi.rel[-D], pseudo.fit, pseudo.refit, thres)
    
    for(k in 1:K){
      
      result.lg.k = result.lg[[k]]
      
      S.list[, -D, k] = result.lg.k[[1]]
      Omega.lg.list[, -D, k] = result.lg.k[[2]]
      Omega.rf.list[, -D, k] = result.lg.k[[3]]
      edge[, -D, k] = result.lg.k[[4]]
      time[k] = result.lg.k[[5]]
      
      record.list[, , k] = result.lg.k[[6]]
      time.list[, k] = result.lg.k[[7]]
    }
    
    result.gg = simu.gglasso(pos, Corr, Sigma, lambda.c, epi.abs[D], epi.rel[D], pseudo.fit, pseudo.refit, thres)
    
    S.list[, D, ] = result.gg[[1]]
    Omega.lg.list[, D, ] = result.gg[[2]]
    Omega.rf.list[, D, ] = result.gg[[3]]
    edge[, D, ] = result.gg[[4]]
  
  }else{
    
    record.list = array(0, c(L, D, K)); time.list = matrix(0, D, K)
    
    result.lg = foreach(k = 1:K, .combine = "list", .multicombine = TRUE, .maxcombine = K, .export = c("simu.lglasso.pos")) %dopar%
      simu.lglasso.pos(pos[k], Corr, Sigma, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres)
    
    for(k in 1:K){
      
      result.lg.k = result.lg[[k]]
      
      S.list[, , k] = result.lg.k[[1]]
      Omega.lg.list[, , k] = result.lg.k[[2]]
      Omega.rf.list[, , k] = result.lg.k[[3]]
      edge[, , k] = result.lg.k[[4]]
      time[k] = result.lg.k[[5]]
      
      record.list[, , k] = result.lg.k[[6]]
      time.list[, k] = result.lg.k[[7]]
    }
  }
  
  return(list(S.list, Omega.lg.list, Omega.rf.list, edge, time, record.list, time.list))  
}


#cross validation selection###################################################################################################################################################################################################################################

#Input###
#fold: number of folds in cross validation
#corr.ind: 0: use Sigma in model fitting, 1: use Corr in model fitting

#Output###
#result: a list of length=fold containing the results from simu.lglasso
#cv: L by D by K by fold by 2 cv scores (the last dimension corresponds to whether we use the adjusted bandwidth h or not)

cv.select = function(pos, h, fold, X, corr.ind, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres){
  
  p = dim(X)[1]; N = dim(X)[2]; K = length(pos); D = length(d.l); L = length(lambda.c)
  
  result = vector("list", fold); cv = array(0, c(L, D, K, fold, 2))
  
  for(i in 1:fold){
    
    pos.test = seq(i, N, fold); pos.train = (1:N)[-pos.test]
    
    Sigma.train = gene.Sigma(X, pos.train, h)
    if(corr.ind == 1){
      Corr.train = gene.corr(X, pos.train, h)
    }else{
      Corr.train = Sigma.train
    }
    Sigma.test = gene.Sigma(X, pos.test, h)
    Sigma.test2 = gene.Sigma(X, pos.test, h*(fold-1)^(1/5))
    
    result.i = simu.lglasso(pos, Corr.train, Sigma.train, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres)
    Omega.rf.list = result.i[[3]]
    result[[i]] = result.i
    
    for(d in 1:D){for(l in L:1){for(k in 1:K){
      
      Omega.rf = as.matrix(Omega.rf.list[[l, d, k]])
      cv[l, d, k, i, 1] = sum(c(t(Sigma.test[, , pos[k]]))*c(Omega.rf)) - log(det(Omega.rf))
      cv[l, d, k, i, 2] = sum(c(t(Sigma.test2[, , pos[k]]))*c(Omega.rf)) - log(det(Omega.rf))
    }}}
    
  }
  
  return(list(result, cv))
}