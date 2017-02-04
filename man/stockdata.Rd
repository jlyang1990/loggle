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

date.index <- rownames(stockdata$data)
stock.sector <- stockdata$info[, "sector"]

# select stocks from sector Materials,
# Telecommunications Services and Utilities
sp.m <- t(stockdata$data[date.index < "2011-01-01", 
stock.sector == "Materials"])
sp.t <- t(stockdata$data[date.index < "2011-01-01", 
stock.sector == "Telecommunications Services"])
sp.u <- t(stockdata$data[date.index < "2011-01-01", 
stock.sector == "Utilities"])
sp <- rbind(sp.m, sp.t, sp.u)

# construct data matrix by taking log ratio of prices
# between adjacent time points
p <- dim(sp)[1]
N <- dim(sp)[2]-1
X <- matrix(0, p, N)
for(i in 1:p) {
  X[i, ] <- scale(log(sp[i, -1] / sp[i, -(N+1)]))
}
dim(X)  # dimension of data matrix

# positions of time points to estimate graphs
pos.example <- round(seq(0.005, 0.995, length=25)
*(ncol(X)-1)+1, 0)
# estimate time-varying graphs and conduct model 
# selection via cross-validation
result <- LGGM.cv(X, pos = pos.example, h = 0.1, 
d.list = c(0, 0.01, 0.05, 0.15, 0.25, 0.35, 1), 
lambda.list = 10 ^ seq(-0.8, -0.2, length = 4), 
fit.type = "pseudo", num.thread = 2)

# graphs at some positions of time points
# in selected model
pos.loc <- round(seq(2, length(pos.example)-1, length = 4))
for(k in 1:length(pos.loc)) {
  adj.matrix <- result$cv.select.result$Omega.edge.list.min[[pos.loc[k]]]
  net <- graph.adjacency(adj.matrix, mode = "undirected", diag = FALSE)
  set.seed(0)
  plot(net, edge.color = "black", vertex.size = 3, vertex.label = NA, 
  vertex.color = rep(rainbow(3), times = c(nrow(sp.m), nrow(sp.t), nrow(sp.u))), 
  layout = layout.circle, main = date.index[pos.example[pos.loc[k]]])
  legend("bottomright", c("Materials", "Tele. Services", "Utilities"), 
  pch = 19, col = rainbow(3), cex = 0.5)
}
}