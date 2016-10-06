\name{LGGM.cv.select}
\alias{LGGM.cv.select}
\title{A function to learn time-varying graphical models with structure smoothness based on cross validation}

\description{
This function is to efficiently conduct model selection based on cross validation for learning time-varying graphical models through a local group-Lasso type penalty.
}

\usage{
LGGM.cv.select <- function(cv.result, 
select.type = "all_flexible", cv.vote.thres = 0.8)
}

\arguments{
  \item{cv.result}{}
  \item{select.type}{}
  \item{cv.vote.thres}{}
}
  
\details{

}

\value{
  \item{d.min}
  \item{lambda.min}
  \item{cv.score.min}
  \item{cv.score.min.sd}
  \item{edge.num.list.min}
  \item{edge.list.min}
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
# estimate time-varying graphs and conduct model selection via cross-validation
result <- LGGM.cv(X, pos = pos.example, fit.type = "pseudo", refit.type = "pseudo", 
h = 0.25, d.list = c(0, 0.01, 0.05, 0.15, 0.25, 0.35, 1), 
lambda.list = c(0.15, 0.2, 0.25, 0.3), num.thread = 2)

# numbers of edges at time point positions in selected model
print(result$cv.select.result$edge.num.list.min)
}