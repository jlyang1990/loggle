
source("TVGM-FUN.R")
load("Data-50.RData") ## X: data matrix; Omega.t: true precision matrices; alpha: generation parameter 


########################
#Data generation#############################################################################################
if(0){
set.seed(1)

p = 50; N = 2001; alpha = 0.28
list[Omega.t, X, edge.t] = gene.data(p, N, alpha)

index = seq(1, 2001, 2)
N = 1001; Omega.t = Omega.t[, , index]; X = X[, index]; edge.t = edge.t[index]
}

####################
#Sigma generation: sufficient statistics ############################################################################################

h = 0.2 ##smoothing parameter
Sigma = gene.Sigma(X, 1:N, h) ## smoothed "sample covariance"
Corr = gene.corr(X, 1:N, h) ## convert Sigma to correlation matrix 

#########################
#Algorithm Parameter initialization####################################################################################

K = 51; pos = (1:K-1)/(K-1)*(N-1)+1  ## positions to fit the model  
d.l = c(0.0001, 0.001, 0.01, seq(0.025, 0.4, 0.025), 1) ## sequence of the neighborhood width parameter d 
lambda.c = seq(0.1, 0.3, length = 21) ## sequence of tuning parameter lambda 
epi.abs = c(rep(1e-5, 11), rep(1e-4, 9)); epi.rel = c(rep(1e-3, 11), rep(1e-2, 9)) ## ADMM stopping rule: for small d , use more strigenet stopping rule



thres = p   ## early stop the grid search of lambda (from largest to smallest) when the fitted model has number of edges > thre*p; when thre=p/2., this early stopping rule is not in effect. 

pseudo.fit = 2 ## 0: glasso; 1: pseudo-likelihood without symmetry; 2: pseudo-likelihood with symmetry; 3: space
pseudo.refit = 1 ## 0: glasso; 1: regression based re-fit


#################################
#Model fitting ############################################################################################

registerDoParallel(4) ##for each : number of threads 

time.ini = proc.time()

list[S.list, Omega.lg.list, Omega.rf.list, edge, time, record.list, time.list] = simu.lglasso(pos, Corr, Sigma, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres)

time.overall = proc.time() - time.ini

stopImplicitCluster()

#####################################
##evaluation FDR, power : 21by 20 by 2, [,,2] are standard deviations over K positions 
list[overlap.edge, FDR, power] = FP(pos, Omega.rf.list, edge, Omega.t, edge.t)

save(S.list, Omega.lg.list, Omega.rf.list, FDR, power, edge, overlap.edge, time, time.overall, record.list, time.list, file = "TVGM-SIMU-50.RData")


#CV Simulation Study#########################################################################################

registerDoParallel(16)

list[result, cv] = cv.select(pos, h, fold = 5, X, corr.ind = 0, d.l, lambda.c, epi.abs, epi.rel, pseudo.fit, pseudo.refit, thres)

stopImplicitCluster()

save(result, cv, file = "TVGM-SIMU-CV-50.RData")