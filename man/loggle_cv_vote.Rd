\name{loggle.cv.vote}
\alias{loggle.cv.vote}
\title{A function to learn time-varying graphical models using cv.vote}

\description{
This function is to efficiently learn time-varying graphical models when hyperparameters are pre-specified. Compared with \link{loggle}, the method cv.vote is implemented. The coefficients in time-varying graphs can be estimated via model refitting.
}

\usage{
loggle.cv.vote(X, pos = 1:ncol(X), 
h = 0.8*ncol(X)^(-1/5), d = 0.2, lambda = 0.25, 
cv.fold = 5, fit.type = c("likelihood", "pseudo", 
"space"), refit = TRUE, cv.vote.thres = 0.8, 
epi.abs = 1e-5, epi.rel = 1e-3, max.step = 500, 
detrend = TRUE, fit.corr = TRUE, num.thread = 1, 
print.detail = TRUE)
}

\arguments{
  \item{X}{a p by N data matrix containing observations on a time grid: p -- number of variables, N -- number of time points}
  \item{pos}{a vector containing elements from 1, 2, ..., N: positions of time points where graphs are estimated}
  \item{h}{a scalar between 0 and 1: bandwidth in kernel estimated sample covariance/correlation matrix (assume time points lie in [0, 1]), default = 0.8*ncol(X)^(-1/5)}
  \item{d}{a scalar or a vector of the same length as \code{pos} with values between 0 and 1: width of neighborhood centered at each position specified by \code{pos}. When it is a scalar, it is shared by each position, default = 0.2}
  \item{lambda}{a scalar or a vector of the same length as \code{pos}: tuning parameter of lasso penalty at each position specified by \code{pos}. When it is a scalar, it is shared by each position, default = 0.25}
  \item{cv.fold}{a scalar: number of cross-validation folds, default = 5}
  \item{fit.type}{a string: "likelihood" -- likelihood estimation, "pseudo" -- pseudo likelihood estimation, or "space" -- sparse partial correlation estimation, default = "pseudo"}
  \item{refit}{logic: if TRUE, conduct additional model refitting after learning graph structures, default = TRUE}
  \item{cv.vote.thres}{a scalar between 0 and 1: an edge is kept after cv.vote if and only if it exists in no less than \code{cv.vote.thres}*\code{cv.fold} cv folds, default = 0.8}
  \item{epi.abs}{a scalar: absolute tolerance in ADMM stopping criterion, default = 1e-5}
  \item{epi.rel}{a scalar: relative tolerance in ADMM stopping criterion, default = 1e-3}
  \item{max.step}{an integer: maximum steps in ADMM iteration, default = 500}
  \item{detrend}{logic: if TRUE, subtract kernel weighted moving average for each variable in data matrix (i.e., detrending), if FALSE, subtract overall average for each variable in data matrix (i.e., centralization), default = TRUE}
  \item{fit.corr}{logic: if TRUE, use sample correlation matrix in model fitting, if FALSE, use sample covariance matrix in model fitting, default = TRUE}
  \item{num.thread}{an integer: number of threads used in parallel computing, default = 1}
  \item{print.detail}{logic: if TRUE, print details in model fitting procedure, default = TRUE}
}
  
\details{
This function calls \link{loggle} in each cross-validation step, and uses cv.vote to get the estimated time-varying graphs from the cross-validation result.

The model fitting method based on pseudo-likelihood (\code{fit.type = "pseudo"} or \code{fit.type = "space"}) is usually less computationally intensive than that based on likelihood (\code{fit.type = "likelihood"}), with similar model fitting performance.

\code{cv.vote.thres} controls the tradeoff between false discovery rate and power. A large value of \code{cv.vote.thres} would decrease false discovery rate but also hurt power.

If no pre-processing has been done to the data matrix \code{X}, \code{detrend = TRUE} is suggested to detrend each variable in data matrix by subtracting corresponding kernel weighted moving average.

\code{fit.corr = TRUE} is suggested in model fitting such that lasso-type penalty can perform better when all the variables are of similar scales. If \code{fit.corr = FALSE} is used, the default value of \code{lambda} should be changed accordingly.
}

\value{
  \item{result.fold}{a list of model fitting results from \link{loggle} for each cv fold}
  \item{Omega}{if refit = TRUE: a list of estimated precision matrices at positions of time points specified by \code{pos}; if refit = FALSE: a list of estimated adjacency matrices at positions of time points specified by \code{pos}}
  \item{edge.num}{a vector of numbers of graph edges at positions of time points specified by \code{pos}}
  \item{edge}{a list of graph edges at positions of time points specified by \code{pos}}
}

\references{
Peng, J., Wang, P., Zhou, N., & Zhu, J. (2012). Partial correlation estimation by joint sparse regression models. Journal of the American Statistical Association.
}

\author{
Yang, J. and Peng, J.
}

\seealso{
\link{loggle} for learning time-varying graphical models, \link{loggle.cv} for learning time-varying graphical models via cross validation, \link{loggle.cv.h} for learning time-varying graphical models via cross validation (with \code{h} fixed), \link{loggle.cv.select} for model selection based on cross validation results, \link{loggle.cv.select.h} for model selection based on cross validation results (with \code{h} fixed), \link{loggle.refit} for model refitting based on estimated graphs.
}

\examples{
data(example)  # load data matrix
dim(X)  # dimension of data matrix
p <- nrow(X)  # number of variables

# positions of time points to estimate graphs
pos <- round(seq(0.02, 0.98, length=49)*(ncol(X)-1)+1)
K <- length(pos)
# estimate time-varying graphs
# num.thread can be set as large as number of cores 
# on a multi-core machine
ts <- proc.time()
result <- loggle.cv.vote(X, pos, h = 0.1, d = 0.15, 
lambda = 0.225, fit.type = "pseudo", refit = TRUE, 
cv.vote.thres = 0.8, num.thread = 1)
te <- proc.time()
sprintf("Time used for loggle.cv.vote: \%.2fs", (te-ts)[3])

# numbers of edges at each time point
print(cbind("time" = seq(0.02, 0.98, length=49),
"edge.num" = result$edge.num))

# graphs at certain time points
par(mfrow = c(2, 4))
pos.plot <- pos[round(seq(1, K, length=16))]
for(k in 1:length(pos.plot)) {
  adj.matrix <- result$Omega[[which(pos == 
  pos.plot[k])]] != 0
  net <- graph.adjacency(adj.matrix, mode = 
  "undirected", diag = FALSE)
  set.seed(0)
  plot(net, vertex.size = 10, vertex.color = 
  "lightblue", vertex.label = NA, edge.color = 
  "black", layout = layout.circle)
  title(main = paste("t =", 
  round(pos.plot[k]/(ncol(X)-1), 2)), cex.main = 0.8)
}

# false discovery rate (FDR) and power based on 
# true precision matrices
edge.num.true <- sapply(1:K, function(i) 
(sum(Omega.true[[pos[i]]]!=0)-p)/2)
edge.num.overlap <- sapply(1:K, function(i) 
(sum(result$Omega[[i]] & Omega.true[[pos[i]]])-p)/2)
perform.matrix <- cbind(
"FDR" = 1 - edge.num.overlap / result$edge.num,
"power" = edge.num.overlap / edge.num.true)
print(apply(perform.matrix, 2, mean))
}