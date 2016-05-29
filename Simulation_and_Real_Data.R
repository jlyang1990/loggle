
source("TVGM-FUN-EVAL.R")

########################
#Simulation study#############################################################################################################################################################################################################################################

########################
#data loading##################################################################################################################

load("Data-500.RData")
load("TVGM-SIMU-500.RData")
load("TVGM-SIMU-500-CV.RData")

K = 25; pos = round(seq(0.02, 0.98, length=K)*(N-1)+1, 0)
d.l = c(0.0001, 0.001, 0.01, seq(0.025, 0.1, 0.025), seq(0.15, 0.3, 0.05), 0.4, 1)
lambda.c = seq(0.17, 0.35, length = 10)

########################
#initial graph check###########################################################################################################

#plot of edge numbers###
plot(seq(0, 1, length=N), edge.t, type = "l", xlab = "Time point", ylab = "Edge number", main = "Plot of edge numbers")

#proportion of overlapping edges with a certain graph###
prop.overlap = rep(0, N-1)
k.list = c(301, 601, 1001)
for(k in k.list){
  index = (1:N)[-k]
  for(i in 1:(N-1)){prop.overlap[i] = (sum(Omega.t[, , k]&Omega.t[, , index[i]])-p)/(sum(Omega.t[, , k]|Omega.t[, , index[i]])-p)}
  plot(seq(0, 1, length=N)[-k], prop.overlap, type = "l", xlab = "Time point", ylab = "Proportion", main = paste("Proportion of overlapping edges with the ", k, "st graph", sep = ""))
  abline(v=(k-1)/(N-1), lty=2)
}

########################
#method comparison#############################################################################################################

#plot of ROC curves###
ROC.plot = function(FDR, power, d.l, index){
  D = dim(FDR)[2]; l = length(index)
  plot(NULL, xlim = c(0, 0.7), ylim = c(0.3, 1), xlab = "FDR", ylab = "power")
  sapply(1:l, function(j) lines(FDR[, index[j], 1], power[, index[j], 1], lwd = 2, col = (rainbow(l))[j]))
  lines(FDR[, 1, 1], power[, 1, 1], lty = 5, lwd = 2, col = "black")
  lines(FDR[, D, 1], power[, D, 1], lty = 4, lwd = 2, col = "grey")
  legend("bottomright", c("d = 0 (Zhou's)", sapply(1:(l-1), function(j) paste("d =", as.character(d.l[index[j]])))), lty = c(5, rep(1, (l-1))), lwd = 2, col = c("black", rainbow(l)[1:(l-1)]), cex = 0.7)
  legend("bottomleft", c(sapply((l-0):l, function(j) paste("d =", as.character(d.l[index[j]]))), "d = 1 (Wang's)"), lty = c(rep(1, 1), 4), lwd = 2, col = c(rainbow(l)[(l-0):l], "grey"), cex = 0.7)
}

ROC.plot(FDR, power, d.l, c(3,5,7,8,9))

#optimal model selection###
#optimal model for our method
list[lambda_d.opt, d.opt, lambda_d.avg.opt, FDR.opt, power.opt, F1.opt] = optimal.select(pos, d.l, lambda.c, d.pos = 1:length(d.l), edge, overlap.edge, edge.t)
#optimal model for Zhou's method
list[lambda_d.opt, d.opt, lambda_d.avg.opt, FDR.opt, power.opt, F1.opt] = optimal.select(pos, d.l, lambda.c, d.pos = 1, edge, overlap.edge, edge.t)
#optimal model for Wang's method
list[lambda_d.opt, d.opt, lambda_d.avg.opt, FDR.opt, power.opt, F1.opt] = optimal.select(pos, d.l, lambda.c, d.pos = length(d.l), edge, overlap.edge, edge.t)

########################
#model selection###############################################################################################################

#model selection via cross validation###
cv[is.na(cv)] = Inf
#model selection for our method
list[S.cv, S.d.cv, S.avg.cv, lambda_d.cv, d.cv, lambda_d.avg.cv, FDR.cv, power.cv, FDR.cv.vote, power.cv.vote, F1.cv, cv.score] = cv.vote(pos, cv[,,,,2], result, d.l, lambda.c, d.pos = 1:length(d.l), Omega.t, edge.t, edge, overlap.edge, vote.thres = 1)
#model selection for Zhou's method
list[S.cv, S.d.cv, S.avg.cv, lambda_d.cv, d.cv, lambda_d.avg.cv, FDR.cv, power.cv, FDR.cv.vote, power.cv.vote, F1.cv, cv.score] = cv.vote(pos, cv[,,,,2], result, d.l, lambda.c, d.pos = 1, Omega.t, edge.t, edge, overlap.edge, vote.thres = 1)
#model selection for Wang's method
list[S.cv, S.d.cv, S.avg.cv, lambda_d.cv, d.cv, lambda_d.avg.cv, FDR.cv, power.cv, FDR.cv.vote, power.cv.vote, F1.cv, cv.score] = cv.vote(pos, cv[,,,,2], result, d.l, lambda.c, d.pos = length(d.l), Omega.t, edge.t, edge, overlap.edge, vote.thres = 1)

#coarse grid search in model selection (large p)###
#data loading
load("TVGM-SIMU-500-CV-s.RData")
K = 9; pos = round(seq(0.02, 0.98, length=K)*(N-1)+1, 0)
d.l = c(0.0001, 0.01, 0.05, 0.15, 0.25, 1)
lambda.c = seq(0.17, 0.29, length = 4)
cv[is.na(cv)] = Inf
#coarse grid search in model selection for our method
list[lambda_d.cv, d.cv, lambda_d.avg.cv, cv.score] = cv.eval(cv[,,,,2], result, d.l, lambda.c, d.pos = 1:length(d.l))
#coarse grid search in model selection for Zhou's method
list[lambda_d.avg.cv, FDR.cv, power.cv, cv.score] = cv.eval(cv[,,,,2], result, d.l, lambda.c, d.pos = 1)
#coarse grid search in model selection for Wang's method
list[lambda_d.cv, d.cv, lambda_d.avg.cv, cv.score] = cv.eval(cv[,,,,2], result, d.l, lambda.c, d.pos = length(d.l))


########################
#Real data application########################################################################################################################################################################################################################################

########################
#data generation###############################################################################################################

library(huge)
data(stockdata)
sp1 = t(stockdata[[1]][, stockdata[[2]][,2]=="Health Care"])
sp2 = t(stockdata[[1]][, stockdata[[2]][,2]=="Information Technology"])
sp3 = t(stockdata[[1]][, stockdata[[2]][,2]=="Telecommunications Services"])
sp4 = t(stockdata[[1]][, stockdata[[2]][,2]=="Utilities"])
sp5 = t(stockdata[[1]][, stockdata[[2]][,2]=="Consumer Discretionary"])
sp6 = t(stockdata[[1]][, stockdata[[2]][,2]=="Consumer Staples"])
sp7 = t(stockdata[[1]][, stockdata[[2]][,2]=="Energy"])
sp8 = t(stockdata[[1]][, stockdata[[2]][,2]=="Financials"])
sp9 = t(stockdata[[1]][, stockdata[[2]][,2]=="Industrials"])
sp10 = t(stockdata[[1]][, stockdata[[2]][,2]=="Materials"])
sp = rbind(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9, sp10)
p = dim(sp)[1]; N = dim(sp)[2]-1
X = matrix(0, p, N)
for(i in 1:p){X[i, ] = scale(log(sp[i, -1]/sp[i, -(N+1)]))}

load("TVGM-SIMU-REAL-CV.RData")

K = 201; pos = round(seq(0.005, 0.995, length=K)*(N-1)+1, 0)
d.l = c(0.0001, 1)
lambda.c = seq(0.3, 0.65, length = 8)
sp.p = c(nrow(sp1), nrow(sp2), nrow(sp3), nrow(sp4), nrow(sp5), nrow(sp6), nrow(sp7), nrow(sp8), nrow(sp9), nrow(sp10))
sp.ind = c(0, cumsum(sp.p))

########################
#initial check#################################################################################################################

#ACF plots for randomly selected stocks to check time independency###
acf(X[1,], main = "ACF plot for 1st stock")
acf(X[91,], main = "ACF plot for 91st stock")
acf(X[181,], main = "ACF plot for 181st stock")
acf(X[272,], main = "ACF plot for 272nd stock")
acf(X[362,], main = "ACF plot for 362nd stock")
acf(X[452,], main = "ACF plot for 452nd stock")

########################
#graph plots for telecom. services and utilities###############################################################################

k.list = c(1, 68, 134, 201)
date.list = c("01/01/2003", "08/09/2004", "04/24/2006", "01/01/2008")

for(k in 1:length(k.list)){
  net = graph.adjacency(Omega.cv[(cumsum(sp.p)[2]+1):cumsum(sp.p)[4], (cumsum(sp.p)[2]+1):cumsum(sp.p)[4], k.list[k]], mode = "undirected", weighted = NULL, diag = FALSE)
  index = rep(NA, 38); index[1:6] = c("AMT", "AT&T", "CTL", "FTR", "Sprint", "Verizon")
  E(net)$color = "black"
  plot(net, vertex.size = 12, vertex.color = c(rep("red",6), rep("blue",32)), layout = layout.fruchterman.reingold, vertex.label = index, main = date.list[k])
  legend("bottomright", c("Telecom. Services", "Utilities"), pch = 20, col = c("red", "blue"), cex = 0.6)
}

########################
#model selection and method comparison#########################################################################################

cv[is.nan(cv)]=Inf

#model selection for Zhou's method
list[S.cv, S.d.cv, S.avg.cv, Omega.cv, Omega.d.cv, Omega.avg.cv, lambda_d.cv, d.cv, lambda_d.avg.cv, cv.score] = cv.real(cv[,,,,2], result, p, N, d.l, lambda.c, d.pos = 1, vote.thres = 0.8)

edge.sector.zhou = matrix(0, 12, K)
for(i in 1:10){edge.sector.zhou[i, ] = sapply(1:K, function(k) (sum(Omega.d.cv[(sp.ind[i]+1):sp.ind[i+1],(sp.ind[i]+1):sp.ind[i+1],k])-sp.p[i])/2)}
edge.sector.zhou[12, ] = sapply(1:K, function(k) (sum(Omega.d.cv[,,k])-p)/2)
edge.sector.zhou[11, ] = edge.sector.zhou[12, ] - colSums(edge.sector.zhou[1:10, ])

#proportion of edges (within separate sectors and between sectors) for Zhou's method
plot(seq(0.005, 0.995, length=K), edge.sector.zhou[11, ]/(p*(p-1)/2-sum(sp.p*(sp.p-1)/2)), type = "l", col = 1, ylim = c(0, 0.5), xlab = "Time", ylab = "Proportion", main = "Proportion of Edges", xaxt = "n")
for(i in 1:10){lines(seq(0.005, 0.995, length=K), edge.sector.zhou[i, ]/(sp.p[i]*(sp.p[i]-1)/2), col = rainbow(10)[i])}
legend("topright", c("Health Care", "I.T.", "Telecom. Services", "Utilities", "Consumer Discretionary", "Consumer Staples", "Energy", "Financials", "Industrials", "Materials", "Between Sectors"), lty = 1, col = c(rainbow(10), 1), cex = 0.55)
axis(side=1, seq(0.005, 0.995, length=4), labels = c("01/01/2003", "08/09/2004", "04/24/2006", "01/01/2008"))

#proportion of edges (within sectors and between sectors) for Zhou's method
plot(seq(0.005, 0.995, length=K), edge.sector.zhou[11, ]/(p*(p-1)/2-sum(sp.p*(sp.p-1)/2)), type = "l", col = 1, ylim = c(0, 0.1), xlab = "Time", ylab = "Proportion", main = "Proportion of Edges", xaxt = "n")
lines(seq(0.005, 0.995, length=K), colSums(edge.sector.zhou[1:10, ])/(sum(sp.p*(sp.p-1)/2)), col = 2)
legend("topright", c("Within Sectors", "Between Sectors"), lty = 1, col = c(2, 1))
axis(side=1, seq(0.005, 0.995, length=4), labels = c("01/01/2003", "08/09/2004", "04/24/2006", "01/01/2008"))

#model selection for Wang's method
list[S.cv, S.d.cv, S.avg.cv, Omega.cv, Omega.d.cv, Omega.avg.cv, lambda_d.cv, d.cv, lambda_d.avg.cv, cv.score] = cv.real(cv[,,,,2], result, p, N, d.l, lambda.c, d.pos = 2, vote.thres = 0.8)

edge.sector.wang = matrix(0, 12, K)
for(i in 1:10){edge.sector.wang[i, ] = sapply(1:K, function(k) (sum(Omega.d.cv[(sp.ind[i]+1):sp.ind[i+1],(sp.ind[i]+1):sp.ind[i+1],k])-sp.p[i])/2)}
edge.sector.wang[12, ] = sapply(1:K, function(k) (sum(Omega.d.cv[,,k])-p)/2)
edge.sector.wang[11, ] = edge.sector.wang[12, ] - colSums(edge.sector.wang[1:10, ])

#proportion of edges (within separate sectors and between sectors) for Wang's method
plot(seq(0.005, 0.995, length=K), edge.sector.wang[11, ]/(p*(p-1)/2-sum(sp.p*(sp.p-1)/2)), type = "l", col = 1, ylim = c(0, 0.5), xlab = "Time", ylab = "Proportion", main = "Proportion of Edges", xaxt = "n")
for(i in 1:10){lines(seq(0.005, 0.995, length=K), edge.sector.wang[i, ]/(sp.p[i]*(sp.p[i]-1)/2), col = rainbow(10)[i])}
legend("topright", c("Health Care", "I.T.", "Telecom. Services", "Utilities", "Consumer Discretionary", "Consumer Staples", "Energy", "Financials", "Industrials", "Materials", "Between Sectors"), lty = 1, col = c(rainbow(10), 1), cex = 0.55)
axis(side=1, seq(0.005, 0.995, length=4), labels = c("01/01/2003", "08/09/2004", "04/24/2006", "01/01/2008"))

#proportion of edges (within sectors and between sectors) for Wang's method
plot(seq(0.005, 0.995, length=K), edge.sector.wang[11, ]/(p*(p-1)/2-sum(sp.p*(sp.p-1)/2)), type = "l", col = 1, ylim = c(0, 0.1), xlab = "Time", ylab = "Proportion", main = "Proportion of Edges", xaxt = "n")
lines(seq(0.005, 0.995, length=K), colSums(edge.sector.wang[1:10, ])/(sum(sp.p*(sp.p-1)/2)), col = 2)
legend("topright", c("Within Sectors", "Between Sectors"), lty = 1, col = c(2, 1))
axis(side=1, seq(0.005, 0.995, length=4), labels = c("01/01/2003", "08/09/2004", "04/24/2006", "01/01/2008"))

#plot of edge numbers for Zhou's and Wang's methods
plot(seq(0.005, 0.995, length=K), edge.sector.wang[12, ], type = "l", ylim = c(0, 1300), xlab = "Time", ylab = "Edge Number", main = "Plot of Edge Numbers", xaxt = "n")
lines(seq(0.005, 0.995, length=K), edge.sector.zhou[12, ], type = "l", lty = 2)
legend("topright", c("Wang's", "Zhou's"), lty = c(1, 2))
axis(side=1, seq(0.005, 0.995, length=4), labels = c("01/01/2003", "08/09/2004", "04/24/2006", "01/01/2008"))

#proportion of edges for Zhou's and Wang's methods (using number of detected edges as denominator)
plot(seq(0.005, 0.995, length=K), colSums(edge.sector.wang[1:10, ])/edge.sector.wang[12, ], type = "l", ylim = c(0, 1.2), xlab = "Time", ylab = "Proportion", main = "Proportion of Edges", xaxt = "n", col = 2)
lines(seq(0.005, 0.995, length=K), edge.sector.wang[11, ]/edge.sector.wang[12, ], type = "l")
lines(seq(0.005, 0.995, length=K), colSums(edge.sector.zhou[1:10, ])/edge.sector.zhou[12, ], type = "l", lty = 2, col = 2)
lines(seq(0.005, 0.995, length=K), edge.sector.zhou[11, ]/edge.sector.zhou[12, ], type = "l", lty = 2)
legend("topright", c("Wang's Within Sectors", "Wang's Between Sectors", "Zhou's Within Sectors", "Zhou's Between Sectors"), lty = c(1, 1, 2, 2), col = c(2, 1, 2, 1), cex = 0.6)
axis(side=1, seq(0.005, 0.995, length=4), labels = c("01/01/2003", "08/09/2004", "04/24/2006", "01/01/2008"))








#parameter estimation###
cv.coef = function(pos, S, Sigma, Omega.t){
  
  K = length(pos); p = dim(Sigma)[1]; N = dim(Sigma)[3]; rho = 0.6; epi.abs = 1e-4; epi.rel = 1e-2
  
  d.KL = rep(0, K); d.entropy = rep(0, K); d.L2 = rep(0, K)
  
  d.table = matrix(0, 3, 2)
  
  for(k in 1:K){
    
    #Omega.rf = refit(pos[k], Sigma, S[[k, drop=F]], matrix(0, p, p), 0.6, epi.abs, epi.rel, m.iter)
    
    Z = rep(0, p*p); U = rep(0, p*p); S.k = c(t(S[[k]])); S.k.L = nrow(S[[k]]); Sigma.k = Sigma[, , pos[k]]
    
    result = .C("ADMM_refit",
                as.integer(p),
                as.integer(0),
                as.double(c(Sigma.k)),
                Z = as.double(Z),
                as.double(U),
                as.integer(S.k),
                as.integer(S.k.L),
                as.double(rho),
                as.double(epi.abs),
                as.double(epi.rel)
    )
    
    Z = matrix(result$Z, p, p)
    Omega.rf = Z+t(Z); diag(Omega.rf) = diag(Omega.rf)/2
    Omega.tk = Omega.t[, , pos[k]]
    Sigma.rf = solve(Omega.rf); Sigma.tk = solve(Omega.tk)
    logdet.Omega.rf = log(det(Omega.rf)); logdet.Omega.tk = log(det(Omega.tk))
    
    d.KL[k] = sum(c(Omega.rf)*c(Sigma.tk)) - logdet.Omega.rf + logdet.Omega.tk - p
    
    d.entropy[k] = sum(c(Sigma.rf)*c(Omega.tk)) + logdet.Omega.rf - logdet.Omega.tk - p
    
    d.L2[k] = sqrt(sum((Sigma.rf - Sigma.tk)^2)/sum(Sigma.tk^2))
  }
  
  d.table[, 1] = c(mean(d.KL), mean(d.entropy), mean(d.L2))
  d.table[, 2] = c(sd(d.KL), sd(d.entropy), sd(d.L2))/sqrt(K)
  
  return(list(d.KL, d.entropy, d.L2, d.table))
}

Sigma = gene.Sigma(X, 1:N, h = 0.2)

S.t = vector("list", K)
for(k in 1:K){S = which(Omega.t[, , pos[k]]!=0, arr.ind=T); S.t[[k]] = S[(S[, 1] - S[, 2])>0, ]}  

#distances for our method
list[d.KL, d.entropy, d.L2, d.table] = cv.coef(pos, S.cv, Sigma, Omega.t)
list[d.KL, d.entropy, d.L2, d.table] = cv.coef(pos, S.d.cv, Sigma, Omega.t)
list[d.KL, d.entropy, d.L2, d.table] = cv.coef(pos, S.avg.cv, Sigma, Omega.t)
#distances for Zhou's method
list[d.KL, d.entropy, d.L2, d.table] = cv.coef(pos, S.d.cv, Sigma, Omega.t)
list[d.KL, d.entropy, d.L2, d.table] = cv.coef(pos, S.avg.cv, Sigma, Omega.t)
#distances for Wang's method
list[d.KL, d.entropy, d.L2, d.table] = cv.coef(pos, S.d.cv, Sigma, Omega.t)
list[d.KL, d.entropy, d.L2, d.table] = cv.coef(pos, S.avg.cv, Sigma, Omega.t)
#distances for gold standard
list[d.KL, d.entropy, d.L2, d.table] = cv.coef(pos, S.t, Sigma, Omega.t)

#distances for sample estimator
d.KL.s = rep(0, K); d.entropy.s = rep(0, K); d.L2.s = rep(0, K)

for(k in 1:K){
  
  Sigma.rf = Sigma[, , pos[k]]; Omega.tk = Omega.t[, , pos[k]]
  Omega.rf = solve(Sigma.rf); Sigma.tk = solve(Omega.tk)
  logdet.Omega.rf = log(det(Omega.rf)); logdet.Omega.tk = log(det(Omega.tk))
  
  d.KL.s[k] = sum(c(Omega.rf)*c(Sigma.tk)) - logdet.Omega.rf + logdet.Omega.tk - p
  
  d.entropy.s[k] = sum(c(Sigma.rf)*c(Omega.tk)) + logdet.Omega.rf - logdet.Omega.tk - p
  
  d.L2.s[k] = sqrt(sum((Sigma.rf - Sigma.tk)^2)/sum(Sigma.tk^2))
}

d.table = matrix(0, 3, 2)
d.table[, 1] = c(mean(d.KL.s), mean(d.entropy.s), mean(d.L2.s))
d.table[, 2] = c(sd(d.KL.s), sd(d.entropy.s), sd(d.L2.s))/sqrt(K)