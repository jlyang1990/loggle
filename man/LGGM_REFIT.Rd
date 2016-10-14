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
pos.example <- round(seq(0.02, 0.98, length=25)*(ncol(X)-1)+1, 0)
# estimate time-varying graphs via cross-validation
result <- LGGM.cv(X, pos = pos.example, fit.type = "pseudo", refit.type = "pseudo", 
h = 0.25, d.list = c(0, 0.01, 0.05, 0.15, 0.25, 0.35, 1), 
lambda.list = c(0.15, 0.2, 0.25, 0.3), return.select = FALSE, num.thread = 2)

# conduct model selection using cross-validation results
cv.select.result <- LGGM.cv.select(cv.result = result, 
select.type = "all_flexible", cv.vote.thres = 0.6)

# adjacency matrices at each position of time point in selected model
Omega.edge.list.min <- cv.select.result$Omega.edge.list.min
# Precision matrices at each position of time point in selected model via model refitting
Omega.rf.list.min <- LGGM.refit(X, pos, Omega.edge.list.min)
}