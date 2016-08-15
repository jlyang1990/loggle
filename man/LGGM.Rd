\name{LGGM}
\alias{LGGM}
\title{A function to learn time-varying graphical models through a local group-Lasso type penalty}
\description{
This function is to efficiently implement the ADMM (alternating directions method of multipliers) algorithm for learning time-varying graphical models through a local group-Lasso type penalty. The objective function to be optimized can be based on negative log-likelihood or pseudo log-likelihood. It can also improve the parameter estimation in time-varying graphical models through model refitting.
}
\usage{
LGGM(X, pos = 1:ncol(X), fit.type = c("glasso", "pseudo", "space"),
refit.type = c("likelihood", "pseudo"), h = 0.8*ncol(X)^(-1/5), d, lambda,
epi.abs = 1e-5, epi.rel = 1e-3, fit.corr = TRUE, num.core = 1)
}

\arguments{
  \item{X}{a p by N data matrix: p -- number of variables, N -- number of time points (sample size)}
  \item{pos}{a vector which is a subset of 1, 2, ..., N: positions of time points where graphs are estimated}
  \item{fit.type}{a string: "glasso" -- graphical Lasso estimation, "pseudo" -- pseudo likelihood estimation, or "space" -- sparse partial correlation estimation, default = "glasso"}
  \item{refit.type}{a string: "likelihood" -- likelihood estimation, or "pseudo" -- pseudo likelihood estimation, default = "likelihood"}
  \item{h}{a scalar: bandwidth in kernel estimated sample covariance matrix or sample correlation matrix}
  \item{d}{a scalar or a vector of the same length as "pos": width of neighborhood centered at each position where graph is estimated. When it is a scalar, widths of neighborhoods are the same across all positions}
  \item{lambda}{a scalar or a vector of the same length as "pos": tuning parameter of Lasso penalty at each position where graph is estimated. When it is a scalar, tuning parameters of Lasso penalty are the same across all positions}
  \item{epi.abs}{a scalar: absolute tolerance in ADMM stopping criterion, default = 1e-5}
  \item{epi.rel}{a scalar: relative tolerance in ADMM stopping criterion, default = 1e-3}
  \item{fit.corr}{logic: whether to use sample correlation matrix rather than sample covariance matrix in model fitting, default = TRUE, corresponding to sample correlation matrix}
  \item{num.core}{an integer: number of cores used in parallel computing, default = 1}
  }
\details{
    score implements the hill climbing algorithm. It can be used to build an ensemble of  DAGs 
    (in form of adjacency matrices) based on bootstrap resamples of the data.
}
\value{
  \item{Omega.list}{a list of estimated precision matrices via model fitting at positions of time points specified by "pos"}
  \item{Omega.rf.list}{a list of estimated precision matrices via model refitting at positions of time points specified by "pos"}
  \item{edge.num.list}{a vector of numbers of graph edges at positions of time points specified by "pos"}
  \item{edge.list}{a list of graph edges at positions of time points specified by "pos"}
}
\references{
Yang, J. and Peng, J. (2016). Estimating Time-Varying Graphical Models through a Local Group-Lasso Type Penalty. arXiv:
Peng, J., Wang, P., Zhou, N., & Zhu, J. (2012). Partial correlation estimation by joint sparse regression models. Journal of the American Statistical Association.
}
\author{
Yang, J. and Peng, J.
}

\examples{
  data(example)
  Y.n=example$Y  # data matrix
  true.dir=example$true.dir  # adjacency matrix of the data generating DAG
  true.ske=example$true.ske  # skeleton graph of the data generating DAG
  
  
  temp=score(Y=Y.n, n.boot=0, score.type="BIC") ## learn DAG using "BIC" score
  adj=temp$adj.matrix
  
  sum(adj)  # number of edges 
  
  
  # Find the total number of skeleton edges and the number of correct skeleton edges 
  # from the adjacency matrix  of an estimated DAG.
  # Then compared with the true.ske.
  
    diag(adj)=0
    tt=adj+t(adj) ##symmetrization
    correct.c=sum((tt>0)&(true.ske>0))/2    
    total.c=sum(tt>0)/2
    total.true=sum(true.ske>0)/2
    c(total.c, correct.c, total.true)

 # Build an ensemble of DAGs based on bootstrap resamples of the data 
  temp.boot=score(Y.n, n.boot=10, score.type="BIC")
  boot.adj=temp.boot$adj.matrix
}