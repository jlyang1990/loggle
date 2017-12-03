\name{loggle.refit}
\alias{loggle.refit}
\title{A function to implement model refitting based on graph structures}

\description{
This function is to efficiently implement model refitting based on time-varying graph structures. Time-varying graph structures can be learned from \link{loggle}, \link{loggle.cv} and \link{loggle.cv.h}.
}

\usage{
loggle.refit(X, pos, adj.mat, h = 0.8*ncol(X)^(-1/5))
}

\arguments{
  \item{X}{a p by N data matrix containing observations on a time grid: p -- number of variables, N -- number of time points}
  \item{pos}{a vector containing elements from 1, 2, ..., N: positions of time points where graphs are estimated}
  \item{adj.mat}{a list of the same length as \code{pos}: adjacency matrices at positions of time points specified by \code{pos}}
  \item{h}{a scalar between 0 and 1: bandwidth in kernel estimated sample covariance/correlation matrix (assume time points lie in [0, 1]), default = 0.8*ncol(X)^(-1/5)}
}

\details{
The model refitting method in this function is based on R package \code{glasso}.
}

\value{
  \item{Omega}{a list of estimated precision matrices via model refitting at positions of time points specified by \code{pos}}
}

\references{
Peng, J., Wang, P., Zhou, N., & Zhu, J. (2012). Partial correlation estimation by joint sparse regression models. Journal of the American Statistical Association.
}

\author{
Yang, J. and Peng, J.
}

\seealso{
\link{loggle} for learning time-varying graphical models, \link{loggle.cv} for learning time-varying graphical models via cross validation, \link{loggle.cv.h} for learning time-varying graphical models via cross validation (with \code{h} fixed), \link{loggle.cv.select} for model selection based on cross validation results, \link{loggle.cv.select.h} for model selection based on cross validation results (with \code{h} fixed), \link{loggle.cv.vote} for learning time-varying graphical models using cv.vote.
}

\examples{
data(example)  # load data matrix
dim(X)  # dimension of data matrix

# positions of time points to estimate graphs
pos <- round(seq(0.02, 0.98, length=25)*(ncol(X)-1)+1)
# estimate time-varying graphs and conduct model 
# selection via cross-validation
result <- loggle.cv.h(X, pos, h = 0.2, 
d.list = c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 1), 
lambda.list = c(0.15, 0.2, 0.25, 0.3), fit.type 
= "pseudo", cv.vote.thres = 1, num.thread = 1)

# adjacency matrices at each time point in selected model
adj.mat.opt <- result$cv.select.result$adj.mat.opt
# precision matrices at each time point in selected model 
# via model refitting
Omega.opt <- loggle.refit(X, pos, adj.mat.opt)
}