\name{stockdata}
\alias{stockdata}
\docType{data}
\title{Stock prices of S&P 500 companies from 2007 to 2017}

\description{
This dataset consists of stock closing prices and company information.
}

\usage{
data(stockdata)
}

\format{
  The format is a list containing two matrices:
  
  1. data - 2518x455: 455 stocks' closing prices for 2518 trading days from 2007 to 2017.
  
  2. info - 455x3:
  The 1st column: the ticker symbol for each company. 
  The 2nd column: the sector for each company.
  The 3rd column: the full name of each company.
}

\details{
This dataset can be used to perform high-dimensional graph estimation to analyze the (time-varying) relationships between S&P 500 companies.
}

\source{
It is publicly available at www.yahoo.com.
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