
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


#optimal model selection (F1 score)###########################################################################################################################################################################################################################

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

#cv.vote######################################################################################################################################################################################################################################################

cv.vote = function(pos, cv, result.cv, d.l, lambda.c, d.pos, Omega.t, edge.t, edge, overlap.edge, vote.thres){
  
  d.l = d.l[d.pos]; edge = edge[, d.pos, , drop=F]; overlap.edge = overlap.edge[, d.pos, , drop=F]
  
  p = dim(Omega.t)[1]; N = dim(Omega.t)[3]; L = dim(cv)[1]; D = length(d.pos); K = dim(cv)[3]; fold = length(result.cv)
  
  edge.cv = rep(0, K); overlap.edge.cv = rep(0, K); edge.d.cv = rep(0, K); overlap.edge.d.cv = rep(0, K); edge.avg.cv = rep(0, K); overlap.edge.avg.cv = rep(0, K)
  s.cv = rep(0, K); s.d.cv = rep(0, K); s.avg.cv = rep(0, K)
  
  Omega.cv = array(0, c(p, p, K, fold)); Omega.d.cv = array(0, c(p, p, K, fold)); Omega.avg.cv = array(0, c(p, p, K, fold))
  
  S.cv = vector("list", K); S.d.cv = vector("list", K); S.avg.cv = vector("list", K)
  
  lambda_d.cv = matrix(0, 2, K); d.cv = rep(0, K+1); lambda_d.avg.cv = rep(0, 2)
  
  FDR.cv = matrix(0, 3, 2); power.cv = matrix(0, 3, 2); FDR.cv.vote = matrix(0, 3, 2); power.cv.vote = matrix(0, 3, 2); F1.cv = matrix(0, 3, 2)
  
  cv = apply(cv[, d.pos, , , drop=F], c(1,2,3), sum)
  
  lambda_d.index = matrix(0, 2, K)
  for(k in 1:K){
    index = which(as.matrix(cv[, , k], L, D) == min(cv[, , k]), arr.ind = T)
    if(nrow(index)>1){index = matrix(index[nrow(index), ], c(1, 2))}
    lambda_d.index[, k] = index
    lambda_d.cv[, k] = c(d.l[index[2]], lambda.c[index[1]])
  }
  
  d.index = which.min(sapply(1:D, function(d) sum(apply(cv[, d, ], 2, min))))
  d.lambda.index = sapply(1:K, function(k) which.min(cv[, d.index, k]))
  d.cv[1] = d.l[d.index]; d.cv[-1] = lambda.c[d.lambda.index]
  
  cv.avg = apply(cv, c(1, 2), mean)
  lambda_d.avg.index = which(cv.avg == min(cv.avg), arr.ind = T)
  lambda_d.avg.cv = c(d.l[lambda_d.avg.index[2]], lambda.c[lambda_d.avg.index[1]])
  
  for(k in 1:K){
    
    edge.cv[k] = edge[lambda_d.index[1, k], lambda_d.index[2, k], k]
    overlap.edge.cv[k] = overlap.edge[lambda_d.index[1, k], lambda_d.index[2, k], k]
    s.cv[k] = cv[lambda_d.index[1, k], lambda_d.index[2, k], k]
    edge.d.cv[k] = edge[d.lambda.index[k], d.index, k]
    overlap.edge.d.cv[k] = overlap.edge[d.lambda.index[k], d.index, k]
    s.d.cv[k] = cv[d.lambda.index[k], d.index, k]
  }
  edge.avg.cv = edge[lambda_d.avg.index[1], lambda_d.avg.index[2], ]
  overlap.edge.avg.cv = overlap.edge[lambda_d.avg.index[1], lambda_d.avg.index[2], ]
  s.avg.cv = cv[lambda_d.avg.index[1], lambda_d.avg.index[2], ]
  
  FDR.list = (1-overlap.edge.cv/edge.cv)
  FDR.d.list = (1-overlap.edge.d.cv/edge.d.cv)
  FDR.avg.list = (1-overlap.edge.avg.cv/edge.avg.cv)
  power.list = (overlap.edge.cv/edge.t[pos])
  power.d.list = (overlap.edge.d.cv/edge.t[pos])
  power.avg.list = (overlap.edge.avg.cv/edge.t[pos])
  
  FDR.cv[, 1] = c(mean(FDR.list), mean(FDR.d.list), mean(FDR.avg.list))
  power.cv[, 1] = c(mean(power.list), mean(power.d.list), mean(power.avg.list))
  FDR.cv[, 2] = c(sd(FDR.list), sd(FDR.d.list), sd(FDR.avg.list))/sqrt(K)
  power.cv[, 2] = c(sd(power.list), sd(power.d.list), sd(power.avg.list))/sqrt(K)
  F1.cv[, 1] = c(mean(2*(1-FDR.list)*power.list/(1-FDR.list+power.list)), mean(2*(1-FDR.d.list)*power.d.list/(1-FDR.d.list+power.d.list)), mean(2*(1-FDR.avg.list)*power.avg.list/(1-FDR.avg.list+power.avg.list)))
  cv.score = c(mean(s.cv), mean(s.d.cv), mean(s.avg.cv))
  
  for(i in 1:fold){
    
    Omega.rf.list = result.cv[[i]][[3]][, d.pos, , drop=F]
    
    for(k in 1:K){
      
      Omega.cv[, , k, i] = as.matrix(Omega.rf.list[[lambda_d.index[1, k], lambda_d.index[2, k], k]])
      Omega.d.cv[, , k, i] = as.matrix(Omega.rf.list[[d.lambda.index[k], d.index, k]])
      Omega.avg.cv[, , k, i] = as.matrix(Omega.rf.list[[lambda_d.avg.index[1], lambda_d.avg.index[2], k]])  
    }
  }
  
  Omega.cv = (apply(Omega.cv!=0, c(1,2,3), sum)>=fold*vote.thres) 
  Omega.d.cv = (apply(Omega.d.cv!=0, c(1,2,3), sum)>=fold*vote.thres)
  Omega.avg.cv = (apply(Omega.avg.cv!=0, c(1,2,3), sum)>=fold*vote.thres)
  
  for(k in 1:K){
    
    S = which(Omega.cv[, , k]!=0, arr.ind=T); S.cv[[k]] = S[(S[, 1] - S[, 2])>0, , drop = F]
    S = which(Omega.d.cv[, , k]!=0, arr.ind=T); S.d.cv[[k]] = S[(S[, 1] - S[, 2])>0, , drop = F]
    S = which(Omega.avg.cv[, , k]!=0, arr.ind=T); S.avg.cv[[k]] = S[(S[, 1] - S[, 2])>0, , drop = F]
    
    edge.cv[k] = (sum(Omega.cv[, , k])-p)/2
    overlap.edge.cv[k] = (sum(Omega.cv[, , k]&Omega.t[, , pos[k]])-p)/2
    edge.d.cv[k] = (sum(Omega.d.cv[, , k])-p)/2
    overlap.edge.d.cv[k] = (sum(Omega.d.cv[, , k]&Omega.t[, , pos[k]])-p)/2
    edge.avg.cv[k] = (sum(Omega.avg.cv[, , k])-p)/2
    overlap.edge.avg.cv[k] = (sum(Omega.avg.cv[, , k]&Omega.t[, , pos[k]])-p)/2
  }
  
  FDR.list = (1-overlap.edge.cv/edge.cv); FDR.list[is.nan(FDR.list)] = 0
  FDR.d.list = (1-overlap.edge.d.cv/edge.d.cv); FDR.d.list[is.nan(FDR.d.list)] = 0
  FDR.avg.list = (1-overlap.edge.avg.cv/edge.avg.cv); FDR.avg.list[is.nan(FDR.avg.list)] = 0
  power.list = (overlap.edge.cv/edge.t[pos])
  power.d.list = (overlap.edge.d.cv/edge.t[pos])
  power.avg.list = (overlap.edge.avg.cv/edge.t[pos])
  
  FDR.cv.vote[, 1] = c(mean(FDR.list), mean(FDR.d.list), mean(FDR.avg.list))
  power.cv.vote[, 1] = c(mean(power.list), mean(power.d.list), mean(power.avg.list))
  FDR.cv.vote[, 2] = c(sd(FDR.list), sd(FDR.d.list), sd(FDR.avg.list))/sqrt(K)
  power.cv.vote[, 2] = c(sd(power.list), sd(power.d.list), sd(power.avg.list))/sqrt(K)
  F1.cv[, 2] = c(mean(2*(1-FDR.list)*power.list/(1-FDR.list+power.list)), mean(2*(1-FDR.d.list)*power.d.list/(1-FDR.d.list+power.d.list)), mean(2*(1-FDR.avg.list)*power.avg.list/(1-FDR.avg.list+power.avg.list)))
  
  return(list(S.cv, S.d.cv, S.avg.cv, lambda_d.cv, d.cv, lambda_d.avg.cv, FDR.cv, power.cv, FDR.cv.vote, power.cv.vote, F1.cv, cv.score))
}

#cv.eval (used in coarse grid search for large p)#############################################################################################################################################################################################################

cv.eval = function(cv, result.cv, d.l, lambda.c, d.pos){
  
  d.l = d.l[d.pos]; L = dim(cv)[1]; D = length(d.pos); K = dim(cv)[3]
  
  s.cv = rep(0, K); s.d.cv = rep(0, K); s.avg.cv = rep(0, K)
  
  lambda_d.cv = matrix(0, 2, K); d.cv = rep(0, K+1); lambda_d.avg.cv = rep(0, 2)
  
  cv = apply(cv[, d.pos, , , drop = F], c(1,2,3), sum)
  
  lambda_d.index = matrix(0, 2, K)
  for(k in 1:K){
    index = which(matrix(cv[, , k], L, D) == min(cv[, , k]), arr.ind = T)
    if(nrow(index)>1){index = matrix(index[nrow(index), ], c(1, 2))}
    lambda_d.index[, k] = index
    lambda_d.cv[, k] = c(d.l[index[2]], lambda.c[index[1]])
  }
  
  d.index = which.min(sapply(1:D, function(d) sum(apply(cv[, d, ], 2, min))))
  d.lambda.index = sapply(1:K, function(k) which.min(cv[, d.index, k]))
  d.cv[1] = d.l[d.index]; d.cv[-1] = lambda.c[d.lambda.index]
  
  cv.avg = apply(cv, c(1, 2), mean)
  lambda_d.avg.index = which(cv.avg == min(cv.avg), arr.ind = T)
  lambda_d.avg.cv = c(d.l[lambda_d.avg.index[2]], lambda.c[lambda_d.avg.index[1]])
  
  for(k in 1:K){
    s.cv[k] = cv[lambda_d.index[1, k], lambda_d.index[2, k], k]
    s.d.cv[k] = cv[d.lambda.index[k], d.index, k]
  }
  s.avg.cv = cv[lambda_d.avg.index[1], lambda_d.avg.index[2], ]
  
  cv.score = c(mean(s.cv), mean(s.d.cv), mean(s.avg.cv))
  
  return(list(lambda_d.cv, d.cv, lambda_d.avg.cv, cv.score))
}

#cv.real######################################################################################################################################################################################################################################################

cv.real = function(cv, result.cv, p, N, d.l, lambda.c, d.pos, vote.thres){
  
  d.l = d.l[d.pos]; L = dim(cv)[1]; D = length(d.pos); K = dim(cv)[3]; fold = length(result.cv)
  
  Omega.cv = array(0, c(p, p, K, fold)); Omega.d.cv = array(0, c(p, p, K, fold)); Omega.avg.cv = array(0, c(p, p, K, fold))
  
  S.cv = vector("list", K); S.d.cv = vector("list", K); S.avg.cv = vector("list", K)
  
  lambda_d.cv = matrix(0, 2, K); d.cv = rep(0, K+1); lambda_d.avg.cv = rep(0, 2)
  
  cv.score = matrix(0, 3, 2)
  
  cv.o = cv[, d.pos, , , drop=F]; cv = apply(cv[, d.pos, , , drop=F], c(1,2,3), sum)
  
  lambda_d.index = matrix(0, 2, K)
  for(k in 1:K){
    index = which(matrix(cv[, , k], L, D) == min(cv[, , k]), arr.ind = T)
    if(nrow(index)>1){index = matrix(index[nrow(index), ], c(1, 2))}
    lambda_d.index[, k] = index
    lambda_d.cv[, k] = c(d.l[index[2]], lambda.c[index[1]])
  }
  
  d.index = which.min(sapply(1:D, function(d) sum(apply(cv[, d, ], 2, min))))
  d.lambda.index = sapply(1:K, function(k) which.min(cv[, d.index, k]))
  d.cv[1] = d.l[d.index]; d.cv[-1] = lambda.c[d.lambda.index]
  
  cv.avg = apply(cv, c(1, 2), mean)
  lambda_d.avg.index = which(cv.avg == min(cv.avg), arr.ind = T)
  lambda_d.avg.cv = c(d.l[lambda_d.avg.index[2]], lambda.c[lambda_d.avg.index[1]])
  
  for(i in 1:fold){
    
    Omega.rf.list = result.cv[[i]][[3]][, d.pos, , drop=F]
    
    for(k in 1:K){
      
      Omega.cv[, , k, i] = as.matrix(Omega.rf.list[[lambda_d.index[1, k], lambda_d.index[2, k], k]])
      Omega.d.cv[, , k, i] = as.matrix(Omega.rf.list[[d.lambda.index[k], d.index, k]])
      Omega.avg.cv[, , k, i] = as.matrix(Omega.rf.list[[lambda_d.avg.index[1], lambda_d.avg.index[2], k]])  
    }
  }
  
  Omega.cv = (apply(Omega.cv!=0, c(1,2,3), sum)>=fold*vote.thres) 
  Omega.d.cv = (apply(Omega.d.cv!=0, c(1,2,3), sum)>=fold*vote.thres)
  Omega.avg.cv = (apply(Omega.avg.cv!=0, c(1,2,3), sum)>=fold*vote.thres)
  
  for(k in 1:K){
    
    S = which(Omega.cv[, , k]!=0, arr.ind=T); S.cv[[k]] = S[(S[, 1] - S[, 2])>0, , drop = F]
    S = which(Omega.d.cv[, , k]!=0, arr.ind=T); S.d.cv[[k]] = S[(S[, 1] - S[, 2])>0, , drop = F]
    S = which(Omega.avg.cv[, , k]!=0, arr.ind=T); S.avg.cv[[k]] = S[(S[, 1] - S[, 2])>0, , drop = F]
  }
  
  cv.temp = sapply(1:fold, function(i) sum(sapply(1:K, function(k) cv.o[lambda_d.index[1, k], lambda_d.index[2, k], k, i])))
  cv.score[1, ] = c(sum(cv.temp), sqrt(fold)*sd(cv.temp))
  cv.temp = sapply(1:fold, function(i) sum(sapply(1:K, function(k) cv.o[d.lambda.index[k], d.index, k, i])))
  cv.score[2, ] = c(sum(cv.temp), sqrt(fold)*sd(cv.temp))
  cv.temp = sapply(1:fold, function(i) sum(sapply(1:K, function(k) cv.o[lambda_d.avg.index[1], lambda_d.avg.index[2], k, i])))
  cv.score[3, ] = c(sum(cv.temp), sqrt(fold)*sd(cv.temp))
  
  return(list(S.cv, S.d.cv, S.avg.cv, Omega.cv, Omega.d.cv, Omega.avg.cv, lambda_d.cv, d.cv, lambda_d.avg.cv, cv.score))
}