\name{LGGM.cv}
\alias{LGGM.cv}
\title{A function to learn time-varying graphical models with structure smoothness based on cross validation}

\description{
This function is to efficiently conduct model selection based on cross validation for learning time-varying graphical models through a local group-Lasso type penalty.
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
  \item{d}{a scalar or a vector of the same length as "pos": width of neighborhood centered at each position where graph is estimated. When it is a scalar, widths of neighborhoods are the same across all positions}
  \item{lambda}{a scalar or a vector of the same length as "pos": tuning parameter of Lasso penalty at each position where graph is estimated. When it is a scalar, tuning parameters of Lasso penalty are the same across all positions}
  \item{epi.abs}{a scalar: absolute tolerance in ADMM stopping criterion, default = 1e-5}
  \item{epi.rel}{a scalar: relative tolerance in ADMM stopping criterion, default = 1e-3}
  \item{fit.corr}{logic: whether to use sample correlation matrix rather than sample covariance matrix in model fitting, default = TRUE, corresponding to sample correlation matrix}
  \item{num.thread}{an integer: number of threads used in parallel computing, default = 1}
}
  
\details{
The model fitting method based on pseudo log-likelihood (\code{fit.type = "pseudo"} or \code{fit.type = "space"}) is usually less computationally intensive than the model fitting method based on negative log-likelihood (\code{fit.type = "glasso"}), with similar (or sometimes even better) model fitting performance. The model refitting method based on pseudo log-likelihood (\code{refit.type = "pseudo"}) is also computationally cheaper than the model refitting method based on negative log-likelihood (\code{refit.type = "likelihood"}).
}

\value{
  \item{Omega.list}{a list of estimated precision matrices via model fitting at positions of time points specified by "pos"}
  \item{Omega.rf.list}{a list of estimated precision matrices via model refitting at positions of time points specified by "pos"}
  \item{edge.num.list}{a vector of numbers of graph edges at positions of time points specified by "pos"}
  \item{edge.list}{a list of graph edges at positions of time points specified by "pos"}
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