
source("TVGM-FUN.R")


########################
#Data generation#############################################################################################

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


########################
#Algorithm parameter initialization####################################################################################

h = 0.15

K = 201; pos = round(seq(0.005, 0.995, length=K)*(N-1)+1, 0)
d.l = c(0.0001, 1)
lambda.c = seq(0.3, 0.65, length = 8)
epi.abs = rep(1e-3, 2); epi.rel = rep(1e-2, 2)

thres = p

pseudo.fit = 2
pseudo.refit = 0


########################
#CV simulation study############################################################################################

registerDoParallel(6)

time.ini = proc.time()

list[result, cv] = cv.select(pos, h, fold = 5, X, corr.ind = 0, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres)

time.overall = proc.time() - time.ini

stopImplicitCluster()

save(result, cv, file = "TVGM-SIMU-REAL-CV.RData")