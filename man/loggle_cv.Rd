\name{loggle.cv}
\alias{loggle.cv}
\title{A function to learn time-varying graphical models via cross validation}

\description{
This function is to efficiently conduct model selection via cross validation for learning time-varying graphical models. It conducts grid search for optimal \code{h} (bandwidth in kernel smoothed sample covariance/correlation matrix), \code{d} (width of neighborhood centered at each position where a graph is estimated) and \code{lambda} (tuning parameter of lasso penalty at each position where a graph is estimated).
}

\usage{
loggle.cv(X, pos = 1:ncol(X), 
h.list = seq(0.1, 0.3, 0.05), d.list = c(0, 0.001, 0.01, 
0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 1), 
lambda.list = seq(0.15, 0.35, 0.02), cv.fold = 5, 
fit.type = "pseudo", return.select = TRUE,
select.type = "all_flexible", cv.vote.thres = 0.8, 
early.stop.thres = 5, epi.abs = 1e-4, epi.rel = 1e-2, 
max.step = 500, detrend = TRUE, fit.corr = TRUE, 
h.correct = TRUE, num.thread = 1, print.detail = TRUE)
}

\arguments{
  \item{X}{a p by N data matrix containing observations on a time grid ranging from 0 to 1: p -- number of variables, N -- number of time points. The nominal time for the kth time point is (k-1)/(N-1)}
  \item{pos}{a vector constitutes a subset of \{1, 2, ..., N\}: indices of time points where graphs are estimated, default = 1:N}
  \item{h.list}{a vector with values between 0 and 1: a grid of bandwidths in kernel smoothed sample covariance/correlation matrix, default = seq(0.1, 0.3, 0.05)}
  \item{d.list}{a vector with values between 0 and 1: a grid of widths of neighborhood centered at each time point specified by \code{pos}, default = c(0, 0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 1)}
  \item{lambda.list}{a vector: a grid of tuning parameters of lasso penalty at each time point specified by \code{pos}, default = seq(0.15, 0.35, 0.02)}
  \item{cv.fold}{a scalar: number of cross-validation folds, default = 5}
  \item{fit.type}{a string: "likelihood" -- likelihood estimation, "pseudo" -- pseudo likelihood estimation, or "space" -- sparse partial correlation estimation, default = "pseudo"}
  \item{return.select}{logic: if TRUE, return model selection result, default = TRUE}
  \item{select.type}{a string: "all_flexible" -- optimal \code{d} and \code{lambda} can vary across time points specified by \code{pos}, "d_fixed" -- optimal \code{d} is fixed across time points specified by \code{pos} and optimal \code{lambda} can vary across time points specified by \code{pos}, "all_fixed" -- optimal \code{d} and \code{lambda} are fixed across time points specified by \code{pos}, default = "all_flexible"}
  \item{cv.vote.thres}{a scalar between 0 and 1: an edge is kept after cv.vote if and only if it exists in no less than \code{cv.vote.thres}*\code{cv.fold} cv folds, default = 0.8}
  \item{early.stop.thres}{a scalar: stopping criterion in grid search -- grid search stops when the ratio between number of detected edges and number of variables exceeds \code{early.stop.thres}, default = 5}
  \item{epi.abs}{a scalar or a vector of the same length as \code{d.list}: absolute tolerance in ADMM stopping criterion, default = 1e-4.}
  \item{epi.rel}{a scalar or a vector of the same length as \code{d.list}: relative tolerance in ADMM stopping criterion, default = 1e-2.}
  \item{max.step}{an integer: maximum iteration steps in ADMM iteration, default = 500}
  \item{detrend}{logic: if TRUE, subtract kernel weighted moving average for each variable in data matrix (i.e., detrending), if FALSE, subtract overall average for each variable in data matrix (i.e., centering), default = TRUE}
  \item{fit.corr}{logic: if TRUE, use sample correlation matrix in model fitting, if FALSE, use sample covariance matrix in model fitting, default = TRUE}
  \item{h.correct}{logic: if TRUE, apply bandwidth adjustment for validation sets due to sample size difference, default = TRUE}
  \item{num.thread}{an integer: number of threads used in parallel computing, default = 1}
  \item{print.detail}{logic: if TRUE, print details in model fitting procedure, default = TRUE}
}
  
\details{
This function conducts grid search for optimal \code{h}, \code{d} and \code{lambda}. It is a wrapper of \link{loggle.cv.h} and calls \link{loggle.cv.h} for each candidate \code{h}.

The model fitting method based on pseudo-likelihood (\code{fit.type = "pseudo"} or \code{fit.type = "space"}) is usually less computationally intensive than that based on likelihood (\code{fit.type = "likelihood"}), with similar model fitting performance.

\code{select.type = "all_flexible"} is for the situation where we expect both the extent of structure smoothness (controlled by \code{d}) and the extent of graph sparsity (controlled by \code{lambda}) vary across time points. If only the extent of graph sparsity varies across time points, \code{select.type = "d_fixed"} should be used. If both of them are expected to be similar across time points, \code{select.type = "all_fixed"} should be used.

\code{cv.vote.thres} controls the tradeoff between false discovery rate and power in model selection. A large value of \code{cv.vote.thres} would decrease false discovery rate but also hurt power. 

\code{early.stop.thres} is used as an early stopping criterion. With limited sample size, a very complicated model seldom leads to a competitive cv score. Therefore in grid search for \code{lambda}, when the estimated model is already large determined by \code{early.stop.thres}, we stop fitting models for smaller \code{lambda}'s (as they often lead to more complicated models). In this case, the output corresponding to the remaining \code{lambda}'s will be \code{NA}.

If no pre-processing has been done to the data matrix \code{X}, \code{detrend = TRUE} is recommended to detrend each variable in data matrix by subtracting corresponding kernel weighted moving average.

\code{fit.corr = TRUE} is recommended such that all the variables are of the same scale. If \code{fit.corr = FALSE} is used, the default value of \code{lambda} may need to be changed accordingly.

\code{h.correct = TRUE} is recommended in calculating kernel smoothed sample covariance matrix for validation sets. This is because bandwidth \code{h} should reflect the difference in sample sizes between training and validation sets.
}

\value{
  \item{cv.result.h}{a list of model fitting results from \link{loggle.cv.h} for each \code{h}}
  \item{cv.select.result}{results from \link{loggle.cv.select} if \code{return.select = TRUE}}
}

\references{
Yang, J. & Peng, J. (2018), 'Estimating Time-Varying Graphical Models', arXiv preprint arXiv:1804.03811
}

\author{
Yang, J. and Peng, J.
}

\seealso{
\link{loggle} for learning time-varying graphical models, \link{loggle.cv.select} for model selection based on cross validation results, \link{loggle.cv.vote} for learning time-varying graphical models using cv.vote, \link{loggle.cv.h} for learning time-varying graphical models via cross validation (with \code{h} fixed), \link{loggle.cv.select.h} for model selection based on cross validation results (with \code{h} fixed), \link{loggle.refit} for model refitting based on estimated graphs, \link{example} for description of an example dataset of time-varying graphs, \link{example_source} for generating the example dataset of time-varying graphs, \link{stockdata} for description of stock price dataset of S&P 500 companies, \link{stockdata_source} for obtaining stock price dataset of S&P 500 companies.
}

\examples{
data(example)  # load example dataset
# data matrix and true precision matrices
X <- example$X
Omega.true <- example$Omega.true
dim(X)  # dimension of data matrix
p <- nrow(X)  # number of variables

# positions of time points to estimate graphs
pos <- round(seq(0.02, 0.98, length=25)*(ncol(X)-1)+1)
K <- length(pos)
# estimate time-varying graphs and conduct model 
# selection via cross-validation
# num.thread can be set as large as number of cores 
# on a multi-core machine (however when p is large, 
# memory overflow should also be taken caution of)
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