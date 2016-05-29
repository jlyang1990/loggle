
source("TVGM-FUN.R")
load("Data-500.RData")


########################
#Sigma generation############################################################################################

h = 0.3
Sigma = gene.Sigma(X, 1:N, h)
Corr = gene.corr(X, 1:N, h)


########################
#Algorithm parameter initialization####################################################################################

K = 25; pos = round(seq(0.02, 0.98, length=K)*(N-1)+1, 0)
d.l = c(0.0001, 0.001, 0.01, seq(0.025, 0.1, 0.025), seq(0.15, 0.3, 0.05), 0.4, 1)
lambda.c = seq(0.17, 0.35, length = 10)
epi.abs = rep(1e-3, 13); epi.rel = rep(1e-2, 13)

thres = 1.25

pseudo.fit = 2
pseudo.refit = 1


########################
#Model fitting############################################################################################

registerDoParallel(7)

time.ini = proc.time()

list[S.list, Omega.lg.list, Omega.rf.list, edge, time, record.list, time.list] = simu.lglasso(pos, Corr, Sigma, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres)

time.overall = proc.time() - time.ini

stopImplicitCluster()

list[overlap.edge, FDR, power] = FP(pos, Omega.rf.list, edge, Omega.t, edge.t)

save(S.list, Omega.lg.list, Omega.rf.list, FDR, power, edge, overlap.edge, time, time.overall, record.list, time.list, file = "TVGM-SIMU-500.RData")


########################
#CV simulation Study#########################################################################################

registerDoParallel(7)

list[result, cv] = cv.select(pos, h, fold = 5, X, corr.ind = 1, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres)

stopImplicitCluster()

save(result, cv, file = "TVGM-SIMU-CV-500.RData")


########################
#CV simulation Study (coarse grid search)#####################################################################

K = 9; pos = round(seq(0.02, 0.98, length=K)*(N-1)+1, 0)
d.l = c(0.0001, 0.01, 0.05, 0.15, 0.25, 1)
lambda.c = seq(0.17, 0.29, length = 4)
epi.abs = rep(1e-3, 6); epi.rel = rep(1e-2, 6)

thres = 0.75

pseudo.fit = 2
pseudo.refit = 1

registerDoParallel(7)

list[result, cv] = cv.select(pos, h, fold = 5, X, corr.ind = 1, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres)

stopImplicitCluster()

save(result, cv, file = "TVGM-SIMU-CV-500-s.RData")