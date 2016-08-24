\name{LGGM.cv}
\alias{LGGM.cv}
\title{A function to learn time-varying graphical models with structure smoothness based on cross validation}

\description{
This function is to efficiently conduct model selection based on cross validation for learning time-varying graphical models through a local group-Lasso type penalty. In model selection, the bandwidth in kernel estimated sample covariance/correlation matrix is fixed.
}

\usage{
LGGM.cv <- function(X, pos = 1:ncol(X), fit.type = "glasso", 
refit.type = "likelihood", h = 0.8*ncol(X)^(-1/5), 
d.list = c(0, 0.001, 0.01, 0.025, 0.05, 0.075, 1, 0.15, 0.2, 0.25, 0.3, 0.4, 1), 
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
  \item{num.fold}
  \item{cv.thres}
  \item{return.select}
  \item{select.type}
  \item{cv.vote.thres}
  \item{epi.abs}{a scalar or a vector of the same length as \code{d.list}: absolute tolerance in ADMM stopping criterion, default = 1e-4 if p >= 400 or 1e-5 otherwise. When it is a scalar, absolute tolerances are the same for all elements in \code{d.list}}
  \item{epi.rel}{a scalar or a vector of the same length as \code{d.list}: relative tolerance in ADMM stopping criterion, default = 1e-2 if p >= 400 or 1e-3 otherwise. When it is a scalar, relative tolerances are the same for all elements in \code{d.list}}
  \item{fit.corr}{logic: whether to use sample correlation matrix rather than sample covariance matrix in model fitting, default = TRUE, corresponding to sample correlation matrix}
  \item{h.correct}
  \item{num.thread}{an integer: number of threads used in parallel computing, default = 1}
}
  
\details{
This function conducts grid search for optimal \code{d} and \code{lambda}, while \code{h} is fixed. To figure out the optimal \code{h} as well, one can apply this function to different \code{h}'s or use the function \link{LGGM.cv.h}.
}

\value{
  \item{cv.score}
  \item{cv.result.list}
  \item{cv.select.result}
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
# estimate time-varying graphs
result <- LGGM(X, pos = pos.example, fit.type = "pseudo", 
refit.type = "pseudo", h = 0.25, d = 0.25, lambda = 0.25, num.thread = 7)

print(result$edge.num.list)  # numbers of edges at time point positions