\name{LGGM.refit}
\alias{LGGM.refit}
\title{A function to implement model refitting based on graphical structures}

\description{
This function is to implement model refitting based on graphical structures for learning time-varying graphical models.
}

\usage{
LGGM.refit(X, pos, Omega.edge.list, h = 0.8*ncol(X)^(-1/5))
}

\arguments{
  \item{X}{a p by N data matrix: p -- number of variables, N -- number of time points (sample size)}
  \item{pos}{a vector which is a subset of 1, 2, ..., N: positions of time points where graphs are estimated}
  \item{Omega.edge.list}{a list of the same length as \code{pos}: adjacency matrices at positions of time points specified by \code{pos}}
  \item{h}{a scalar: bandwidth in kernel estimated sample covariance matrix}
}

\details{
This function uses the model refitting method based on negative log-likelihood.
}

\value{
  \item{Omega.rf.list}{a list of estimated precision matrices via model refitting at positions of time points specified by \code{pos}}
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
pos.example <- round(seq(0.02, 0.98, length=25)
*(ncol(X)-1)+1, 0)
# estimate time-varying graphs and conduct model 
# selection via cross-validation
result <- LGGM.cv(X, pos = pos.example, h = 0.2, 
d.list = c(0, 0.01, 0.05, 0.15, 0.25, 0.35, 1), 
lambda.list = c(0.15, 0.2, 0.25, 0.3), fit.type 
= "pseudo", refit.type = "glasso", cv.vote.thres = 1,
num.thread = 2)

# adjacency matrices at each position of time point in 
# selected model
Omega.edge.list.min <- result$cv.select.result
$Omega.edge.list.min
# Precision matrices at each position of time point in 
# selected model via model refitting
Omega.rf.list.min <- LGGM.refit(X, pos.example, 
Omega.edge.list.min)
}