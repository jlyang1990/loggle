\name{LGGM}
\alias{LGGM}
\title{A function to learn time-varying graphical models through a local group-Lasso type penalty}
\description{
This function is to efficiently implement the hill climbing algorithm 
for learning a DAG by optimizing a score.
It can also build an emsemble of DAGs based on bootstrap resamples of the data. 
}
\usage{
score(Y, n.boot=0, score.type="BIC", threshold=0, max.step=500,  ini.adj.matrix=NULL, 
blacklist=NULL, whitelist=NULL, standardize=TRUE,  standardize.boot=TRUE, 
random.forest=FALSE, random.step.length=NULL, 
nrestart=0, perturb=0, shuffle=FALSE, print=FALSE, EPS=1e-06)
}

\arguments{
  \item{Y}{an n by p data matrix: n -- sample size, p -- number of variables}
  \item{n.boot}{an integer: the number of bootstrap resamples of the data matrix Y, 
  default = 0, meaning no bootstrapping}
  \item{score.type}{a string: "BIC" or "likelihood", default = "BIC"} 
 \item{threshold}{a nonnegative scalar: the cutoff value for the change of the score 
 to decide whether to stop the search; 
 default = 0, meaning stop search when score is not improved}
  \item{max.step}{an integer: the maximum number of search steps  of the hill climbing algorithm, default = 500}
  \item{ini.adj.matrix}{ a p by p 0-1 matrix: the initial graph, default = NULL, meaning the empty graph}
  \item{blacklist}{a p by p 0-1 matrix: if the (i,j)th-entry is "1", 
  then the edge i-->j will be excluded from the DAG during the search, default= NULL}
  \item{whitelist}{a p by p 0-1 matrix: if the (i,j)th-entry is "1", 
  then the edge i-->j will always be included in the DAG during the search, default=NULL}
  \item{standardize}{logic: whether to standardize the data to have mean zero and sd one, default = TRUE}
  \item{standardize.boot}{logic: whether to standardize the bootstrap resamples, default = TRUE}
  \item{random.forest}{logic: whether to use the "random forest" idea for further variance reduction, default=FALSE}
  \item{random.step.length}{a vector: specify ``random forest" steps}
  \item{nrestart}{an integer: number of times to restart the search algorithm after a local optimal is achieved. 
   The purpose is to search for global optimal, default=0, meaning no restart.}
  \item{perturb}{an integer: how many random addition/deletion/reversal operations 
  should be used in each random restart,  
  default = 0, corresponding to no restart.}
  \item{shuffle}{logic: whether to shuffle the order of variables before DAG learning. 
  The purpose is  to avoid potential systematic biases in simulation studies 
  (due to  possible coincidence of the topological ordering and the nominal ordering of the variables), default=FALSE}
  \item{print}{logic: whether print the step information, default= FALSE}
  \item{EPS}{a scalar: a number to indicate what small values will be  treated as zero, 
  default = 1e-06, meaning that values with magnitude smaller than 1e-6 will be treated as zero.}
}
\details{
    score implements the hill climbing algorithm. It can be used to build an ensemble of  DAGs 
    (in form of adjacency matrices) based on bootstrap resamples of the data.
}
\value{
    If n.boot=0 (no bootstrap), then a list with five components: 
    \item{adj.matrix}{adjacency matrix for the estimated DAG}
    \item{final.step}{a number recording how many search steps are conducted before the procedure stops}
    \item{movement}{a matrix recording the selected operation, \
    addition, deletion or reversal of an edge, at each search step} 
    \item{delta.min}{a vector recording the change of  score due to the selected operation at each search step} 
    \item{bic.ok}{a vector recording whether there is error in score calculation at each step} 
    If n.boot>0 (with bootstrap) then a list with one component:
    \item{adj.matrix}{adjacency matrices for the estimated DAGs based on n.boot bootstrap  resamples}
}
\references{
Yang, J. and Peng, J. (2016) Estimating Time-Varying Graphical Models through a Local Group-Lasso Type Penalty. arXiv:
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