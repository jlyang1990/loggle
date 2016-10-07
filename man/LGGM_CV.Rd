\name{LGGM.cv}
\alias{LGGM.cv}
\title{A function to learn time-varying graphical models with structure smoothness via cross validation}

\description{
This function is to efficiently conduct model selection via cross validation for learning time-varying graphical models through a local group-Lasso type penalty. In model selection, the bandwidth in kernel estimated sample covariance/correlation matrix should be pre-specified.
}

\usage{
LGGM.cv <- function(X, pos = 1:ncol(X), fit.type = c("glasso", "pseudo", "space"), 
refit.type = c("likelihood", "pseudo"), h = 0.8*ncol(X)^(-1/5), 
d.list = c(0, 0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 1), 
lambda.list = seq(0.15, 0.35, length = 11), num.fold = 5, cv.thres = 5, 
return.select = TRUE, select.type = "all_flexible", cv.vote.thres = 0.8, 
epi.abs = ifelse(nrow(X) >= 400, 1e-4, 1e-5), 
epi.rel = ifelse(nrow(X) >= 400, 1e-2, 1e-3), 
fit.corr = TRUE, h.correct = TRUE, num.thread = 1)
}

\arguments{
  \item{X}{a p by N data matrix: p -- number of variables, N -- number of time points (sample size)}
  \item{pos}{a vector which is a subset of 1, 2, ..., N: positions of time points where graphs are estimated}
  \item{fit.type}{a string: "glasso" -- graphical Lasso estimation, "pseudo" -- pseudo likelihood estimation, or "space" -- sparse partial correlation estimation, default = "glasso"}
  \item{refit.type}{a string: "likelihood" -- likelihood estimation, or "pseudo" -- pseudo likelihood estimation, default = "likelihood"}
  \item{h}{a scalar: bandwidth in kernel estimated sample covariance matrix or sample correlation matrix}
  \item{d.list}{a vector: a grid of widths of neighborhood centered at each position where graph is estimated}
  \item{lambda.list}{a vector: a grid of tuning parameters of Lasso penalty at each position where graph is estimated}
  \item{num.fold}{a scalar: number of cv folds}
  \item{cv.thres}{a scalar: stopping criterion in grid search: grid search stops when the ratio between number of detected edges and number of nodes exceeds \code{cv.thres}}
  \item{return.select}{logic: whether to return model selection result, default = TRUE}
  \item{select.type}{a string: "all_flexible" -- optimal d and lambda can vary at each position where graph is estimated, "d_fixed" -- optimal d is fixed and optimal lambda can vary at each position where graph is estimated, "all_fixed" -- optimal d and lambda are fixed at each position where graph is estimated, default = "all_flexible"}
  \item{cv.vote.thres}{a scalar (between 0 and 1): an edge is retained in cv vote if and only if it exists in no less than \code{cv.vote.thres}*\code{num.fold} cv folds}
  \item{epi.abs}{a scalar or a vector of the same length as \code{d.list}: absolute tolerance in ADMM stopping criterion, default = 1e-4 if p >= 400 or 1e-5 otherwise. When it is a scalar, absolute tolerances are the same for all elements in \code{d.list}}
  \item{epi.rel}{a scalar or a vector of the same length as \code{d.list}: relative tolerance in ADMM stopping criterion, default = 1e-2 if p >= 400 or 1e-3 otherwise. When it is a scalar, relative tolerances are the same for all elements in \code{d.list}}
  \item{fit.corr}{logic: whether to use sample correlation matrix rather than sample covariance matrix in model fitting, default = TRUE, which is sample correlation matrix}
  \item{h.correct}{logic: whether to apply bandwidth correction to h in calculating cv scores for validation sets}
  \item{num.thread}{an integer: number of threads used in parallel computing, default = 1}
}
  
\details{
This function conducts grid search for optimal \code{d} and \code{lambda} selection, while \code{h} is fixed and pre-specified. To figure out the optimal \code{h} as well, one can apply this function to different candidate \code{h}'s or use the function \link{LGGM.cv.h}.

When underlying graphs are sparse, there is usually no need to conduct the complete grid search, especially for those dense graphs which may need a large amount of training time. \code{cv.thres} is used as an early stopping criterion in grid search, where the grid search (grid of \code{lambda} for each \code{d}) stops when the number of detected edges is large enough.

\code{select.type = "all_flexible"} is chosen in model selection if we believe both the extent of smoothness (controlled by \code{d}) and sparsity (controlled by \code{lambda}) vary across time points. If only the extent of sparsity varies across time points, we choose \code{select.type = "d_fix"}. If both of them are consistent across time points, we choose \code{select.type = "all_fix"}.

\code{cv.vote.thres} controls the tradeoff between false discovery rate and power in model selection. A large value of \code{cv.vote.thres} would decrease false discovery rate but also hurt power. 

\code{h.correct = TRUE} is suggested in calculating cross-validation scores for validation sets. The bandwidth correction can make the choice of bandwidth \code{h} consistent between training and validation sets, since the optimal bandwidth \code{h} is related to sample size \code{N}, which differs in training and validation sets. 
}

\value{
  \item{cv.score}{an array of cv scores for each combination of d, lambda, position of time point and cv fold}
  \item{cv.result.list}{a list of model fitting results (of the same format as the results from the function \link{LGGM}) for each cv fold}
  \item{cv.select.result}{results from the function \link{LGGM.cv.select} if \code{return.select = TRUE}}
}

\references{
Peng, J., Wang, P., Zhou, N., & Zhu, J. (2012). Partial correlation estimation by joint sparse regression models. Journal of the American Statistical Association.
}

\author{
Yang, J. and Peng, J.
}

\examples{
data(example)  # load data matrix
dim(X)  # dimension of data matrix

# positions of time points to estimate graphs
pos.example <- round(seq(0.02, 0.98, length=25)*(ncol(X)-1)+1, 0)
# estimate time-varying graphs and conduct model selection via cross-validation
result <- LGGM.cv(X, pos = pos.example, fit.type = "pseudo", refit.type = "pseudo", 
h = 0.25, d.list = c(0, 0.01, 0.05, 0.15, 0.25, 0.35, 1), 
lambda.list = c(0.15, 0.2, 0.25, 0.3), num.thread = 2)

# optimal values of d at each position of time point
print(result$cv.select.result$d.min)
# optimal values of lambda at each position of time point
print(result$cv.select.result$lambda.min)
# numbers of edges at each position of time point in selected model
print(result$cv.select.result$edge.num.list.min)
}