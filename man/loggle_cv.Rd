\name{loggle.cv}
\alias{loggle.cv}
\title{A function to learn time-varying graphical models via cross validation}

\description{
This function is to efficiently conduct model selection via cross validation for learning time-varying graphical models. It conducts grid search for optimal \code{h} (bandwidth in kernel estimated sample covariance/correlation matrix), \code{d} (width of neighborhood centered at each position where graph is estimated) and \code{lambda} (tuning parameters of lasso penalty at each position where graph is estimated) selection.
}

\usage{
loggle.cv <- function(X, pos = 1:ncol(X), 
h.list = seq(0.1, 0.3, 0.05), d.list = c(0, 0.001, 0.01, 
0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 1), 
lambda.list = seq(0.15, 0.35, 0.02), cv.fold = 5, 
fit.type = c("likelihood", "pseudo", "space"), 
return.select = TRUE,
select.type = c("all_flexible", "d_fixed", "all_fixed"), 
cv.vote.thres = 0.8, early.stop.thres = 5, 
epi.abs = 1e-4, epi.rel = 1e-2, max.step = 500, 
detrend = TRUE, fit.corr = TRUE, h.correct = TRUE, 
num.thread = 1, print.detail = TRUE)
}

\arguments{
  \item{X}{a p by N data matrix containing observations on a time grid: p -- number of variables, N -- number of time points}
  \item{pos}{a vector containing elements from 1, 2, ..., N: positions of time points where graphs are estimated}
  \item{h.list}{a vector with values between 0 and 1: a grid of bandwidths in kernel estimated sample covariance/correlation matrix (assume time points lie in [0, 1]), default = seq(0.1, 0.3, 0.05)}
  \item{d.list}{a vector with values between 0 and 1: a grid of widths of neighborhood centered at each position specified by \code{pos}, default = c(0, 0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 1)}
  \item{lambda.list}{a vector: a grid of tuning parameters of lasso penalty at each position specified by \code{pos}, default = seq(0.15, 0.35, 0.02)}
  \item{cv.fold}{a scalar: number of cross-validation folds, default = 5}
  \item{fit.type}{a string: "likelihood" -- likelihood estimation, "pseudo" -- pseudo likelihood estimation, or "space" -- sparse partial correlation estimation, default = "pseudo"}
  \item{return.select}{logic: if TRUE, return model selection result, default = TRUE}
  \item{select.type}{a string: "all_flexible" -- optimal \code{d} and \code{lambda} can vary across positions specified by \code{pos}, "d_fixed" -- optimal \code{d} is fixed and optimal \code{lambda} can vary across positions specified by \code{pos}, "all_fixed" -- optimal \code{d} and \code{lambda} are fixed across positions specified by \code{pos}, default = "all_flexible"}
  \item{cv.vote.thres}{a scalar between 0 and 1: an edge is kept after cv.vote if and only if it exists in no less than \code{cv.vote.thres}*\code{cv.fold} cv folds, default = 0.8}
  \item{early.stop.thres}{a scalar: stopping criterion in grid search -- grid search stops when the ratio between number of detected edges and number of nodes exceeds \code{early.stop.thres}, default = 5}
  \item{epi.abs}{a scalar or a vector of the same length as \code{d.list}: absolute tolerance in ADMM stopping criterion, default = 1e-4. When it is a scalar, it is shared by each element in \code{d.list}}
  \item{epi.rel}{a scalar or a vector of the same length as \code{d.list}: relative tolerance in ADMM stopping criterion, default = 1e-2. When it is a scalar, it is shared by each element in \code{d.list}}
  \item{max.step}{an integer: maximum steps in ADMM iteration, default = 500}
  \item{detrend}{logic: if TRUE, subtract kernel weighted moving average for each variable in data matrix (i.e., detrending), if FALSE, subtract overall average for each variable in data matrix (i.e., centralization), default = TRUE}
  \item{fit.corr}{logic: if TRUE, use sample correlation matrix in model fitting, if FALSE, use sample covariance matrix in model fitting, default = TRUE}
  \item{h.correct}{logic: if TRUE, apply bandwidth correction to \code{h} in calculating cv scores for validation sets, default = TRUE}
  \item{num.thread}{an integer: number of threads used in parallel computing, default = 1}
  \item{print.detail}{logic: if TRUE, print details in model fitting procedure, default = TRUE}
}
  
\details{
This function conducts grid search for optimal \code{h}, \code{d} and \code{lambda} selection. It is a wrapper of \link{loggle.cv.h} and calls \link{loggle.cv.h} for each candidate \code{h}.

The model fitting method based on pseudo-likelihood (\code{fit.type = "pseudo"} or \code{fit.type = "space"}) is usually less computationally intensive than that based on likelihood (\code{fit.type = "likelihood"}), with similar model fitting performance.

\code{select.type = "all_flexible"} is chosen in model selection if we believe both the extent of structure smoothness (controlled by \code{d}) and sparsity (controlled by \code{lambda}) vary across time points. If only the extent of sparsity varies across time points, we choose \code{select.type = "d_fixed"}. If both of them are the same across time points, we choose \code{select.type = "all_fixed"}.

\code{cv.vote.thres} controls the tradeoff between false discovery rate and power in model selection. A large value of \code{cv.vote.thres} would decrease false discovery rate but also hurt power. 

When underlying graphs are sparse, there is usually no need to conduct the complete grid search, especially for those dense graphs which may need a large amount of training time. \code{early.stop.thres} is used as an early stopping criterion in grid search, where the grid search of \code{lambda} for each \code{d} stops when the number of detected edges is large enough. In this case, the output corresponding to the remaining \code{lambda}'s will be \code{NA}.

If no pre-processing has been done to the data matrix \code{X}, \code{detrend = TRUE} is suggested to detrend each variable in data matrix by subtracting corresponding kernel weighted moving average.

\code{fit.corr = TRUE} is suggested in model fitting such that lasso-type penalty can perform better when all the variables are of similar scales. If \code{fit.corr = FALSE} is used, the default value of \code{lambda} should be changed accordingly.

\code{h.correct = TRUE} is suggested in calculating cross-validation scores for validation sets. The bandwidth correction can make the choice of bandwidth \code{h} consistent between training and validation sets, since the optimal bandwidth \code{h} is related to sample size \code{N}, which differs in training and validation sets. 
}

\value{
  \item{cv.result.h}{a list of model fitting results from \link{loggle.cv.h} for each h}
  \item{cv.select.result}{results from \link{loggle.cv.select} if \code{return.select = TRUE}}
}

\references{
Peng, J., Wang, P., Zhou, N., & Zhu, J. (2012). Partial correlation estimation by joint sparse regression models. Journal of the American Statistical Association.
}

\author{
Yang, J. and Peng, J.
}

\seealso{
\link{loggle} for learning time-varying graphical models, \link{loggle.cv.h} for learning time-varying graphical models via cross validation (with \code{h} fixed), \link{loggle.cv.select} for model selection based on cross validation results, \link{loggle.cv.select.h} for model selection based on cross validation results (with \code{h} fixed), \link{loggle.refit} for model refitting based on estimated graphs.
}

\examples{
data(example)  # load data matrix
dim(X)  # dimension of data matrix
p <- nrow(X)  # number of variables

# positions of time points to estimate graphs
pos <- round(seq(0.02, 0.98, length=25)*(ncol(X)-1)+1)
K <- length(pos)
# estimate time-varying graphs and conduct model 
# selection via cross-validation
# num.thread can be set as large as number of cores 
# on a multi-core machine
ts <- proc.time()
result <- loggle.cv(X, pos, 
h.list = c(0.15, 0.2, 0.25, 0.3), 
d.list = c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 1), 
lambda.list = c(0.15, 0.2, 0.25, 0.3), fit.type 
= "pseudo", cv.vote.thres = 1, num.thread = 1)
te <- proc.time()
sprintf("Time used for loggle.cv: \%.2fs", (te-ts)[3])

# optimal values of h, d and lambda, and number of 
# selected edges at each time point
select.result <- result$cv.select.result
print(cbind("time" = seq(0.02, 0.98, length=25),
"h.opt" = rep(select.result$h.opt, K),
"d.opt" = select.result$d.opt,
"lambda.opt" = select.result$lambda.opt,
"edge.num.opt" = select.result$edge.num.opt))

# false discovery rate (FDR) and power based on 
# true precision matrices for selected model
edge.num.opt <- select.result$edge.num.opt
edge.num.true <- sapply(1:K, function(i) 
(sum(Omega.true[[pos[i]]]!=0)-p)/2)
edge.num.overlap <- sapply(1:K, function(i) 
(sum(select.result$adj.mat.opt[[i]]
& Omega.true[[pos[i]]])-p)/2)
perform.matrix <- cbind(
"FDR" = 1 - edge.num.overlap / edge.num.opt,
"power" = edge.num.overlap / edge.num.true)
print(apply(perform.matrix, 2, mean))
}