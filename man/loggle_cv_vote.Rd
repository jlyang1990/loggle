\name{loggle.cv.vote}
\alias{loggle.cv.vote}
\title{A function to learn time-varying graphical models using cv.vote}

\description{
This function is to efficiently learn time-varying graphical models for a given set of tuning parameters. Different from \link{loggle}, cv.vote is also applied in estimating time-varying graphs.
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
  \item{X}{a p by N data matrix containing observations on a time grid ranging from 0 to 1: p -- number of variables, N -- number of time points. The nominal time for the kth time point is (k-1)/(N-1)}
  \item{pos}{a vector constitutes a subset of \{1, 2, ..., N\}: indices of time points where graphs are estimated, default = 1:N}
  \item{h}{a scalar between 0 and 1: bandwidth in kernel smoothed sample covariance/correlation matrix, default = 0.8*N^(-1/5)}
  \item{d}{a scalar or a vector of the same length as \code{pos} with values between 0 and 1: width of neighborhood centered at each time point specified by \code{pos}, default = 0.2}
  \item{lambda}{a scalar or a vector of the same length as \code{pos}: tuning parameter of lasso penalty at each time point specified by \code{pos}, default = 0.25}
  \item{cv.fold}{a scalar: number of cross-validation folds, default = 5}
  \item{fit.type}{a string: "likelihood" -- likelihood estimation, "pseudo" -- pseudo likelihood estimation, or "space" -- sparse partial correlation estimation, default = "pseudo"}
  \item{refit}{logic: if TRUE, conduct model refitting given learned graph structures, default = TRUE}
  \item{cv.vote.thres}{a scalar between 0 and 1: an edge is kept after cv.vote if and only if it exists in no less than \code{cv.vote.thres}*\code{cv.fold} cv folds, default = 0.8}
  \item{epi.abs}{a scalar: absolute tolerance in ADMM stopping criterion, default = 1e-5}
  \item{epi.rel}{a scalar: relative tolerance in ADMM stopping criterion, default = 1e-3}
  \item{max.step}{an integer: maximum steps in ADMM iteration, default = 500}
  \item{detrend}{logic: if TRUE, subtract kernel weighted moving average for each variable in data matrix (i.e., detrending), if FALSE, subtract overall average for each variable in data matrix (i.e., centering), default = TRUE}
  \item{fit.corr}{logic: if TRUE, use sample correlation matrix in model fitting, if FALSE, use sample covariance matrix in model fitting, default = TRUE}
  \item{num.thread}{an integer: number of threads used in parallel computing, default = 1}
  \item{print.detail}{logic: if TRUE, print details in model fitting procedure, default = TRUE}
}
  
\details{
The idea of cross-validation is implemented in this function, where \link{loggle} is applied on each cross-validation fold to get fold-wise estimated time-varying graphs, and cv.vote is then applied across cross-validation folds to get the unified version of estimated time-varying graphs.

The model fitting method based on pseudo-likelihood (\code{fit.type = "pseudo"} or \code{fit.type = "space"}) is usually less computationally intensive than that based on likelihood (\code{fit.type = "likelihood"}), with similar model fitting performance.

\code{cv.vote.thres} controls the tradeoff between false discovery rate and power. A large value of \code{cv.vote.thres} would decrease false discovery rate but also hurt power.

If no pre-processing has been done to the data matrix \code{X}, \code{detrend = TRUE} is recommended to detrend each variable in data matrix by subtracting corresponding kernel weighted moving average.

\code{fit.corr = TRUE} is recommended such that all the variables are of the same scales. If \code{fit.corr = FALSE} is used, the default value of \code{lambda} may need to be changed accordingly.
}

\value{
  \item{result.fold}{a list of model fitting results from \link{loggle} for each cv fold}
  \item{Omega}{a list of estimated precision matrices at time points specified by \code{pos}}
  \item{edge.num}{a vector of numbers of edges at time points specified by \code{pos}}
  \item{edge}{a list of edges at time points specified by \code{pos}}
}

\references{
}

\author{
Yang, J. and Peng, J.
}

\seealso{
\link{loggle} for learning time-varying graphical models, \link{loggle.cv} for learning time-varying graphical models via cross validation, \link{loggle.cv.select} for model selection based on cross validation results.
}

\examples{
data(example)  # load example dataset
# data matrix and true precision matrices
X <- example$X
Omega.true <- example$Omega.true
dim(X)  # dimension of data matrix
p <- nrow(X)  # number of variables

# positions of time points to estimate graphs
pos <- round(seq(0.02, 0.98, length=49)*(ncol(X)-1)+1)
K <- length(pos)
# estimate time-varying graphs
# num.thread can be set as large as number of cores 
# on a multi-core machine (however when p is large, 
# memory overflow should also be taken caution of)
ts <- proc.time()
result <- loggle.cv.vote(X, pos, h = 0.1, d = 0.15, 
lambda = 0.225, fit.type = "pseudo", refit = TRUE, 
cv.vote.thres = 0.8, num.thread = 1)
te <- proc.time()
sprintf("Time used for loggle.cv.vote: \%.2fs", (te-ts)[3])

# number of edges at each time point
print(cbind("time" = seq(0.02, 0.98, length=49),
"edge.num" = result$edge.num))

# graphs at specified time points
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