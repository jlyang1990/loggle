\name{stockdata}
\alias{stockdata}
\docType{data}
\title{Stock price of S&P 500 companies from 2007 to 2017}

\description{
This data set consists of stock price and company information.
}

\usage{data(stockdata)}

\format{
  The format is a list containing conatins two matrices.
  1. data - 1258x452, represents the 452 stocks' close prices for 1258 trading days.
  2. info - 452x3:
  The 1st column: the query symbol for each company. 
  The 2nd column: the categoriy for each company.
  The 3rd column: the full name of each company.
}

\details{
This data set can be used to perform high-dimensional graph estimation to analyze the relationships between S&P 500 companies.
}

\source{
It is publicly available at yahoo.com.
}

\references{
Peng, J., Wang, P., Zhou, N., & Zhu, J. (2012). Partial correlation estimation by joint sparse regression models. Journal of the American Statistical Association.
}

\author{
Yang, J. and Peng, J.
}

\examples{
data(stockdata)
image(stockdata$data)
stockdata$info
}