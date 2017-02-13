\name{example_source}
\alias{example_source}
\docType{data}
\title{Source code for generating underlying true time-varying graphs}

\description{
This is the source code for generating underlying true time-varying graphs (covariance matrices and observed samples).
}

\references{
Peng, J., Wang, P., Zhou, N., & Zhu, J. (2012). Partial correlation estimation by joint sparse regression models. Journal of the American Statistical Association.
}

\author{
Yang, J. and Peng, J.
}

\examples{
library(sparseMVN)
library(matrixcalc)

# function to generate data in time-varying graphs
geneData <- function(p, N, alpha = 0.28) {
  
  done <- FALSE
  
  while(!done) {
    
    # generate base time-varying lower triangular matrices
    t <- seq(0, 1, length = N)
    B1 <- matrix(rnorm(p^2, 0, 1/2), p, p); B1[upper.tri(B1)] <- 0
    B2 <- matrix(rnorm(p^2, 0, 1/2), p, p); B2[upper.tri(B2)] <- 0
    B3 <- matrix(rnorm(p^2, 0, 1/2), p, p); B3[upper.tri(B3)] <- 0
    B4 <- matrix(rnorm(p^2, 0, 1/2), p, p); B4[upper.tri(B4)] <- 0
    uni <- runif(4, -0.5, 0.5)
    
    Omega.true.list <- vector("list", N)
    X <- matrix(0, p, N)
    
    for(i in 1:N) {
      
      # generate raw precision matrix at time point i
      G <- (B1*sin(pi*t[i]/2+uni[1]) + B2*cos(pi*t[i]/2+uni[2]) 
            + B3*sin(pi*t[i]/4+uni[3]) + B4*cos(pi*t[i]/4+uni[4]))/2
      Omega <- G %*% t(G)
      Omega <- Omega %*% diag(1/sqrt((1:p)))
      Omega[upper.tri(Omega)] <- 0
      Omega <- Omega + t(Omega)
      Omega.diag <- diag(Omega)/(2*sqrt(1:p)) + log(p, 10)/4
      
      # implement soft-thresholding to off-diagonals in precision matrix
      Omega.t1 <- matrix(1, p, p) - alpha/abs(Omega)
      Omega.t2 <- matrix(1, p, p) - alpha/(2*abs(Omega))
      Omega <- Omega.t2 * (Omega.t1 > 0) * Omega
      diag(Omega) <- Omega.diag
      
      if(!is.positive.definite(Omega)) {
        break
      }
      
      Omega.true.list[[i]] <- Matrix(Omega, sparse = T)
      
      # generate observed sample from precision matrix
      X[, i] <- rmvn.sparse(1, rep(0, p), Cholesky(Omega.true.list[[i]]))
      
      if(i == N) {
        done = TRUE
      }
    }
  }
  
  result <- new.env()
  result$Omega.true.list <- Omega.true.list
  result$X <- X
  result <- as.list(result)
  
  return(result)
}

# implement geneData function
set.seed(1)
result.data <- geneData(p = 50, N = 2001)
Omega.true.list <- result.data$Omega.true.list[seq(1, 2001, 2)]
X <- result.data$X[, seq(1, 2001, 2)]

# save Omega.true.list and X
save(Omega.true.list, X, file = "example.rda")
}