\name{example}
\alias{example}
\docType{data}
\title{An example dataset of time-varying graphs}

\description{
This dataset contains an example of time-varying graphs including precision matrices and observations, with number of variables being 50.
}

\usage{
data(example)
}

\format{
  This dataset consists of two parts:
  
  1. X - 50x1001: a data matrix containing 1001 observations of 50 variables on an equally-spaced time grid ranging from 0 to 1.
  
  2. Omega.true - a list of length 1001 with each element being a 50x50 matrix: a list of precision matrices on an equally-spaced time grid ranging from 0 to 1.
}

\details{
This dataset can be used for evaluating method performance via comparing the estimated time-varying graphs with the underlying true time-varying graphs. See examples in \link{loggle}, \link{loggle.cv.vote}, \link{loggle.cv}, \link{loggle.cv.h} for more details.
}

\source{
The source code is in \link{example_source}.
}

\references{
}

\author{
Yang, J. and Peng, J.
}