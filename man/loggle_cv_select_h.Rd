\name{loggle.cv.select.h}
\alias{loggle.cv.select.h}
\title{A function to conduct model selection based on cross validation results (with \code{h} fixed)}

\description{
This function is to conduct model selection for time-varying graphical models based on cross validation results from \link{loggle.cv.h}.
}

\usage{
loggle.cv.select.h(cv.result, select.type = "all_flexible", 
cv.vote.thres = 0.8)
}

\arguments{
  \item{cv.result}{a list: results from \link{loggle.cv.h}}
  \item{select.type}{a string: "all_flexible" -- optimal \code{d} and \code{lambda} can vary across time points specified by \code{pos}, "d_fixed" -- optimal \code{d} is fixed across time points specified by \code{pos} and optimal \code{lambda} can vary across time points specified by \code{pos}, "all_fixed" -- optimal \code{d} and \code{lambda} are fixed across time points specified by \code{pos}, default = "all_flexible"}
  \item{cv.vote.thres}{a scalar between 0 and 1: an edge is kept after cv.vote if and only if it exists in no less than \code{cv.vote.thres}*\code{cv.fold} cv folds, default = 0.8}
}
  
\details{
\code{select.type = "all_flexible"} is for the situation where we expect both the extent of structure smoothness (controlled by \code{d}) and the extent of graph sparsity (controlled by \code{lambda}) vary across time points. If only the extent of graph sparsity varies across time points, \code{select.type = "d_fixed"} should be used. If both of them are expected to be similar across time points, \code{select.type = "all_fixed"} should be used.

\code{cv.vote.thres} controls the tradeoff between false discovery rate and power in model selection. A large value of \code{cv.vote.thres} would decrease false discovery rate but also hurt power.
}

\value{
  \item{d.opt}{a vector of optimal values of \code{d} for each estimated graph}
  \item{lambda.opt}{a vector of optimal values of \code{lambda} for each estimated graph}
  \item{cv.score.opt}{optimal cv score (averaged over time points and cv folds)}
  \item{edge.num.opt}{a vector of numbers of edges for each estimated graph}
  \item{edge.opt}{a list of edges for each estimated graph}
  \item{adj.mat.opt}{a list of adjacency matrices for each estimated graph}
}

\references{
Yang, J. & Peng, J. (2018), 'Estimating Time-Varying Graphical Models', arXiv preprint arXiv:1804.03811
}

\author{
Yang, J. and Peng, J.
}

\seealso{
\link{loggle} for learning time-varying graphical models, \link{loggle.cv} for learning time-varying graphical models via cross validation, \link{loggle.cv.select} for model selection based on cross validation results.
}

\examples{
data(example)  # load example dataset
X <- example$X  # data matrix
dim(X)  # dimension of data matrix

# positions of time points to estimate graphs
pos <- round(seq(0.1, 0.9, length=9)*(ncol(X)-1)+1)
# estimate time-varying graphs via cross-validation
result <- loggle.cv.h(X, pos, h = 0.2, 
d.list = c(0, 0.05, 0.15, 1), lambda.list 
= c(0.2, 0.25), cv.fold = 3, fit.type = "pseudo", 
cv.vote.thres = 1, num.thread = 1)

# conduct model selection using cross-validation results
select.result <- loggle.cv.select.h(cv.result = result, 
select.type = "all_flexible", cv.vote.thres = 0.8)

# optimal values of d and lambda, and number of 
# selected edges at each time point
print(cbind("time" = seq(0.1, 0.9, length=9),
"d.opt" = select.result$d.opt,
"lambda.opt" = select.result$lambda.opt,
"edge.num.opt" = select.result$edge.num.opt))
}