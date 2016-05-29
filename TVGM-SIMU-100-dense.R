
source("TVGM-FUN.R")
load("Data-100-dense.RData")


########################
#Data generation#############################################################################################

if(0){
set.seed(1)

p = 100; N = 2001; alpha = 0.28
list[Omega.t, X, edge.t] = gene.data(p, N, alpha)

index = seq(1, 2001, 2)
N = 1001; Omega.t = Omega.t[, , index]; X = X[, index]; edge.t = edge.t[index]
}


########################
#Sigma generation############################################################################################

h = 0.3
Sigma = gene.Sigma(X, 1:N, h)
Corr = gene.corr(X, 1:N, h)


########################
#Algorithm parameter initialization####################################################################################

K = 49; pos = round(seq(0.02, 0.98, length=K)*(N-1)+1, 0)
d.l = c(0.0001, 0.001, 0.01, seq(0.025, 0.1, 0.025), seq(0.15, 0.3, 0.05), 0.4, 1)
lambda.c = seq(0.15, 0.35, length = 11)
epi.abs = c(rep(1e-5, 11), rep(1e-4, 7)); epi.rel = c(rep(1e-3, 11), rep(1e-2, 7))

thres = 3

pseudo.fit = 0
pseudo.refit = 0


########################
#Model fitting############################################################################################

registerDoParallel(7)

time.ini = proc.time()

list[S.list, Omega.lg.list, Omega.rf.list, edge, time, record.list, time.list] = simu.lglasso(pos, Corr, Sigma, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres)

time.overall = proc.time() - time.ini

stopImplicitCluster()

list[overlap.edge, FDR, power] = FP(pos, Omega.rf.list, edge, Omega.t, edge.t)

save(S.list, Omega.lg.list, Omega.rf.list, FDR, power, edge, overlap.edge, time, time.overall, record.list, time.list, file = "TVGM-SIMU-100-dense.RData")


########################
#CV simulation study############################################################################################

registerDoParallel(7)

time.ini = proc.time()

list[result, cv] = cv.select(pos, h, fold = 5, X, corr.ind = 1, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres)

time.overall = proc.time() - time.ini

stopImplicitCluster()

save(result, cv, time.overall, file = "TVGM-SIMU-CV-100-dense.RData")