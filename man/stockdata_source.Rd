\name{stockdata_source}
\alias{stockdata_source}
\docType{data}
\title{Source code for stock prices of S&P 500 companies from 2007 to 2016}

\description{
This is the source code for generating the dataset containing stock closing prices and company information of S&P 500 companies in \link{stockdata}.
}

\details{
Since the components of S&P 500 companies are kept updating, it is recommended that we conduct thorough examination and modification to this source code before we produce the result.
}

\references{
}

\author{
Yang, J. and Peng, J.
}

\examples{
library(XML)
library(RCurl)
library(quantmod)

# extract stock information from wiki
url <- "https://en.wikipedia.org/wiki/List_of_S\%26P_500_companies"
tabs <- getURL(url)
stock.info <- readHTMLTable(tabs, stringsAsFactors = FALSE)[[1]]
stock.ticker <- stock.info$`Ticker symbol`
stock.name <- stock.info$Security
stock.sector <- stock.info$`GICS Sector`

# extract indices of dates and stock closing prices
date.index <- as.character(index(getSymbols(stock.ticker[1], 
auto.assign=FALSE, from='2007-01-01', to="2017-01-01")))
stock.closeprice <- sapply(stock.ticker, function(x) 
as.vector(getSymbols(x, auto.assign=FALSE, from='2007-01-01', 
to="2017-01-01")[, 4]))

# delete stocks not covering the whole time interval
stock.count <- sapply(stock.closeprice, length)
stock.index <- which(stock.count == max(stock.count))
stock.ticker <- stock.ticker[stock.index]
stock.name <- stock.name[stock.index]
stock.sector <- stock.sector[stock.index]
# reformulate matrix of stock closing prices
stock.closeprice <- sapply(stock.closeprice[stock.index], 
function(x) x)
rownames(stock.closeprice) <- date.index

# combine stock information
stock.info <- cbind(stock.ticker, stock.sector, stock.name)
colnames(stock.info) <- c("ticker", "sector", "name")
# delete duplicated stocks
# 22: "GOOGL", 23: "GOOG", 407: "FOXA", 408: "FOX"
stockdata <- list(data = stock.closeprice[, -c(23, 408)], 
info = stock.info[-c(23, 408), ])

# save stockdata into stockdata.rda
save(stockdata, file = "stockdata.rda")
}