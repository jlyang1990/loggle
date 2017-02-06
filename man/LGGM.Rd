\name{LGGM}
\alias{LGGM}
\title{A function to learn time-varying graphical models with structure smoothness}

\description{
This function is to efficiently implement the ADMM (alternating directions method of multipliers) algorithm for learning time-varying graphical models through a local group-Lasso type penalty. The objective function to be optimized is based on negative log-likelihood or pseudo log-likelihood. The parameters in time-varying graphical models can be estimated via model refitting.
}

\usage{
LGGM(X, pos = 1:ncol(X), h = 0.8*ncol(X)^(-1/5), 
d = 0.2, lambda = 0.25, fit.type = c("likelihood", 
"pseudo", "space"), refit = TRUE, epi.abs = 1e-5, 
epi.rel = 1e-3, max.step = 500, detrend = TRUE, 
fit.corr = TRUE, num.thread = 1, print.detail = TRUE)
}

\arguments{
  \item{X}{a p by N data matrix: p -- number of variables, N -- number of time points (sample size)}
  \item{pos}{a vector which is a subset of 1, 2, ..., N: positions of time points where graphs are estimated}
  \item{h}{a scalar between 0 and 1: bandwidth in kernel estimated sample covariance matrix or sample correlation matrix, with the assumption that time points lie in the range from 0 to 1, default = 0.8*ncol(X)^(-1/5)}
  \item{d}{a scalar or a vector of the same length as \code{pos} with values between 0 and 1: width of neighborhood centered at each position where graph is estimated, with the assumption that time points lie in the range from 0 to 1. When it is a scalar, widths of neighborhoods are the same across all positions, default = 0.2}
  \item{lambda}{a scalar or a vector of the same length as \code{pos}: tuning parameter of Lasso penalty at each position where graph is estimated. When it is a scalar, tuning parameters of Lasso penalty are the same across all positions, default = 0.25}
  \item{fit.type}{a string: "likelihood" -- likelihood estimation, "pseudo" -- pseudo likelihood estimation, or "space" -- sparse partial correlation estimation, default = "pseudo"}
  \item{refit}{logic: if TRUE, conduct additional model refitting after graph structures are learned, default = TRUE}
  \item{epi.abs}{a scalar: absolute tolerance in ADMM stopping criterion, default = 1e-5}
  \item{epi.rel}{a scalar: relative tolerance in ADMM stopping criterion, default = 1e-3}
  \item{max.step}{an integer: maximum steps in ADMM iteration, default = 500}
  \item{detrend}{logic: if TRUE, subtract kernel weighted moving average for each variable in data matrix (i.e., detrending), if FALSE, subtract overall average for each variable in data matrix (i.e., centralization), default = TRUE}
  \item{fit.corr}{logic: if TRUE, use sample correlation matrix in model fitting, if FALSE, use sample covariance matrix in model fitting, default = TRUE}
  \item{num.thread}{an integer: number of threads used in parallel computing, default = 1}
  \item{print.detail}{logic: if TRUE, print details in model fitting procedure, default = TRUE}
}
  
\details{
The model fitting method based on pseudo log-likelihood (\code{fit.type = "pseudo"} or \code{fit.type = "space"}) is usually less computationally intensive than the model fitting method based on negative log-likelihood (\code{fit.type = "likelihood"}), with similar (or sometimes even better) model fitting performance.

If no pre-processing has been done to the data matrix \code{X}, \code{detrend = TRUE} is required to detrend each variable in data matrix by subtracting each kernel weighted moving average.

\code{fit.corr = TRUE} is suggested in model fitting since the Lasso-type penalty can be applied appropriately when all the variables are of similar scales. If \code{fit.corr = FALSE} is used, the default value of \code{lambda} should change accordingly.
}

\value{
  \item{Omega.list}{if refit = TRUE: a list of estimated precision matrices via model refitting at positions of time points specified by \code{pos}; if refit = FALSE: a list of estimated precision matrices via model fitting at positions of time points specified by \code{pos}}
  \item{edge.num.list}{a vector of numbers of graph edges at positions of time points specified by \code{pos}}
  \item{edge.list}{a list of graph edges at positions of time points specified by \code{pos}}
}

\references{
Peng, J., Wang, P., Zhou, N., & Zhu, J. (2012). Partial correlation estimation by joint sparse regression models. Journal of the American Statistical Association.
}

\author{
Yang, J. and Peng, J.
}

\seealso{
\link{LGGM.cv} for learning time-varying graphical models via cross validation, \link{LGGM.cv.h} for learning time-varying graphical models via cross validation (with h selection), \link{LGGM.cv.select} for model selection based on cross validation results, and \link{LGGM.refit} for model refitting based on estimated graphs.
}

\examples{
data(example)  # load data matrix
dim(X)  # dimension of data matrix
p <- nrow(X)  # number of variables

# positions of time points to estimate graphs
pos.example <- round(seq(0.02, 0.98, length=49)
*(ncol(X)-1)+1, 0)
K <- length(pos.example)
# estimate time-varying graphs
result <- LGGM(X, pos = pos.example, h = 0.1, 
d = 0.15, lambda = 0.225, fit.type = "pseudo", 
refit = TRUE, num.thread = 2)

# numbers of edges at each position of time point
print(cbind("time" = seq(0.02, 0.98, length=49),
"edge.num" = result$edge.num.list))

# graphs at some positions of time points
par(mfrow = c(2, 4))
pos.plot <- pos.example[round(seq(1, K, length=16))]
for(k in 1:length(pos.plot)) {
  adj.matrix <- result$Omega.list[[which(pos.example 
  == pos.plot[k])]] != 0
  net <- graph.adjacency(adj.matrix, mode = "undirected", 
  diag = FALSE)
  set.seed(0)
  plot(net, vertex.size = 10, vertex.color = "lightblue", 
  vertex.label = NA, edge.color = "black", layout = 
  layout.circle)
  title(main = paste("t =", round(pos.plot[k]/(ncol(X)-1), 
  2)), cex.main = 0.8)
}

# false discovery rate (FDR) and power based on 
# true precision matrices
edge.num.true.list <- sapply(1:K, function(i) 
(sum(Omega.true.list[[pos.example[i]]]!=0)-p)/2)
edge.num.overlap.list <- sapply(1:K, 
function(i) (sum(result$Omega.list[[i]]
& Omega.true.list[[pos.example[i]]])-p)/2)
perform.matrix <- cbind(
"FDR" = 1 - edge.num.overlap.list / result$edge.num.list,
"power" = edge.num.overlap.list / edge.num.true.list)
print(apply(perform.matrix, 2, mean))
}