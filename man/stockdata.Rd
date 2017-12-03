\name{stockdata}
\alias{stockdata}
\docType{data}
\title{Stock prices of S&P 500 companies from 2007 to 2017}

\description{
This dataset consists of stock closing prices and company information of S&P 500 companies.
}

\usage{
data(stockdata)
}

\format{
  This dataset contains two matrices:
  
  1. data - 2518x455: 455 stocks' closing prices for 2518 trading days from 2007 to 2017.
  
  2. info - 455x3:
  1st column: ticker symbol for each company. 
  2nd column: sector for each company.
  3rd column: full name of each company.
}

\details{
This dataset can be used to perform time-varying graph estimations to reveal the dynamic relationships between S&P 500 companies.
}

\source{
Data are publicly available at www.yahoo.com, and R package \code{quantmod} can be used to extract the stock information. The source code can be found in \link{stockdata_source}.
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
pos <- round(seq(0.005, 0.995, length=25)*(N-1)+1)
# estimate time-varying graphs and conduct model 
# selection via cross-validation
result <- loggle.cv.h(X, pos, h = 0.1, 
d.list = c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 1), 
lambda.list = 10 ^ seq(-0.8, -0.2, length = 4), 
fit.type = "pseudo", num.thread = 1)

# graphs at certain time points in selected model
pos.plot <- round(seq(2, length(pos)-1, length = 4))
for(k in 1:length(pos.plot)) {
  adj.matrix <- 
  result$cv.select.result$adj.mat.opt[[pos.plot[k]]]
  net <- graph.adjacency(adj.matrix, mode = "undirected", 
  diag = FALSE)
  set.seed(0)
  V(net)$color <- rep(rainbow(3), times = c(nrow(sp.m), 
  nrow(sp.t), nrow(sp.u)))
  sp.ind <- rep(1:3, times = c(nrow(sp.m), nrow(sp.t), 
  nrow(sp.u)))
  E(net)$color <- apply(get.edgelist(net), 1, function(x) 
  ifelse(sp.ind[x[1]] == sp.ind[x[2]], "black", "gray"))
  E(net)$width <- apply(get.edgelist(net), 1, function(x) 
  ifelse(sp.ind[x[1]] == sp.ind[x[2]], 1, 0.5))
  plot(net, vertex.size = 3, vertex.label = NA, layout = 
  layout.circle, main = date.index[pos[pos.plot[k]]])
  legend("bottomright", c("Materials", "Tele. Services", 
  "Utilities"), pch = 19, col = rainbow(3), cex = 0.4)
  legend("bottomleft", c("within sector", 
  "between sector"), lwd = c(1, 0.4), col = c("black", 
  "gray"), cex = 0.4)
}
}