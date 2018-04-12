\name{loggle.refit}
\alias{loggle.refit}
\title{A function to conduct model refitting given learned graph structures}

\description{
This function is to efficiently conduct model refitting given learned time-varying graph structures. Time-varying graph structures can be learned from \link{loggle}, \link{loggle.cv.vote}, \link{loggle.cv} or \link{loggle.cv.h}.
}

\usage{
loggle.refit(X, pos, adj.mat, h = 0.8*ncol(X)^(-1/5), 
print.detail = TRUE)
}

\arguments{
  \item{X}{a p by N data matrix containing observations on a time grid ranging from 0 to 1: p -- number of variables, N -- number of time points. The nominal time for the kth time point is (k-1)/(N-1)}
  \item{pos}{a vector constitutes a subset of \{1, 2, ..., N\}: indices of time points where graphs are estimated, default = 1:N}
  \item{adj.mat}{a list of the same length as \code{pos}: adjacency matrices at time points specified by \code{pos}}
  \item{h}{a scalar between 0 and 1: bandwidth in kernel smoothed sample covariance/correlation matrix, default = 0.8*N^(-1/5)}
  \item{print.detail}{logic: if TRUE, print details in model refitting procedure, default = TRUE}
}

\details{
Function "glasso" in R package \code{glasso} is applied in model refitting.
}

\value{
  \item{Omega}{a list of estimated precision matrices at time points specified by \code{pos}}
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
pos <- round(seq(0.02, 0.98, length=25)*(ncol(X)-1)+1)
# estimate time-varying graphs and conduct model 
# selection via cross-validation
result <- loggle.cv.h(X, pos, h = 0.2, 
d.list = c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 1), 
lambda.list = c(0.15, 0.2, 0.25, 0.3), fit.type 
= "pseudo", cv.vote.thres = 1, num.thread = 1)

# estimated adjacency matrices at each time point
adj.mat.opt <- result$cv.select.result$adj.mat.opt
# estimated precision matrices at each time point
# via model refitting
Omega.opt <- loggle.refit(X, pos, adj.mat.opt)
}