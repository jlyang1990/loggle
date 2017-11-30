\name{loggle.cv.select.h}
\alias{loggle.cv.select.h}
\title{A function to conduct model selection based on cross validation results (with h selection)}

\description{
This function is to conduct model selection based on cross validation results from function \link{loggle.cv.h} for learning time-varying graphical models.
}

\usage{
loggle.cv.select.h <- function(cv.result, 
select.type = c("all_flexible", "d_fixed", "all_fixed"), 
cv.vote.thres = 0.8)
}

\arguments{
  \item{cv.result}{a list: results from the function \link{loggle.cv.h}}
  \item{select.type}{a string: "all_flexible" -- optimal d and lambda can vary at each position where graph is estimated, "d_fixed" -- optimal d is fixed and optimal lambda can vary at each position where graph is estimated, "all_fixed" -- optimal d and lambda are fixed at each position where graph is estimated, default = "all_flexible"}
  \item{cv.vote.thres}{a scalar (between 0 and 1): an edge is retained in cv vote if and only if it exists in no less than \code{cv.vote.thres}*\code{num.fold} cv folds, default = 0.8}
}
  
\details{
\code{select.type = "all_flexible"} is chosen in model selection if we believe both the extent of smoothness (controlled by \code{d}) and sparsity (controlled by \code{lambda}) vary across time points. If only the extent of sparsity varies across time points, we choose \code{select.type = "d_fixed"}. If both of them are consistent across time points, we choose \code{select.type = "all_fixed"}.

\code{cv.vote.thres} controls the tradeoff between false discovery rate and power in model selection. A large value of \code{cv.vote.thres} would decrease false discovery rate but also hurt power.
}

\value{
  \item{h.min}{optimal value of h}
  \item{d.min}{a vector of optimal values of d at each position of time point}
  \item{lambda.min}{a vector of optimal values of lambda at each position of time point}
  \item{cv.score.min.h}{optimal cv scores for each h}
  \item{cv.score.min}{optimal cv score for optimal h (averaged over positions of time points and cv folds)}
  \item{cv.score.min.sd}{standard deviation of optimal cv scores across cv folds for optimal h}
  \item{edge.num.list.min}{a vector of numbers of graph edges in selected model at each position of time point}
  \item{edge.list.min}{a list of graph edges in selected model at each position of time point}
  \item{Omega.edge.list.min}{a list of adjacency matrices in selected model at each position of time point}
}

\references{
Peng, J., Wang, P., Zhou, N., & Zhu, J. (2012). Partial correlation estimation by joint sparse regression models. Journal of the American Statistical Association.
}

\author{
Yang, J. and Peng, J.
}

\seealso{
\link{loggle} for learning time-varying graphical models, \link{loggle.cv} for learning time-varying graphical models via cross validation, \link{loggle.cv.h} for learning time-varying graphical models via cross validation (with h selection), \link{loggle.cv.select} for model selection based on cross validation results, \link{loggle.refit} for model refitting based on estimated graphs.
}

\examples{
data(example)  # load data matrix
dim(X)  # dimension of data matrix

# positions of time points to estimate graphs
pos.example <- round(seq(0.02, 0.98, length=25)
*(ncol(X)-1)+1, 0)
# list of candidate values of h
h.list <- c(0.15, 0.2, 0.25, 0.3)
# estimate time-varying graphs and conduct model 
# selection via cross-validation
result <- loggle.cv.h(X, pos = pos.example, h.list = h.list, 
d.list = c(0, 0.01, 0.05, 0.15, 0.25, 0.35, 1), 
lambda.list = c(0.15, 0.2, 0.25, 0.3), fit.type = "pseudo", 
cv.vote.thres = 1, num.thread = 1)

# conduct model selection using cross-validation results
cv.select.result <- loggle.cv.select.h(cv.result = result, 
select.type = "all_flexible", cv.vote.thres = 0.8)

# optimal values of d, lambda and number of selected edges
# at each position of time point
print(cbind("time" = seq(0.02, 0.98, length=25),
"h.opt" = rep(cv.select.result$h.min, length(pos.example)),
"d.opt" = cv.select.result$d.min,
"lambda.opt" = cv.select.result$lambda.min,
"edge.num.opt" = cv.select.result$edge.num.list.min))
}