\name{score_shd}
\alias{score_shd}
\title{A function to use structural hamming distance to aggregate DAGs}
\description{This function uses a family of gSHD (generalized structural hamming distance) 
to aggregate an ensemble of DAGs  (e.g., learned on bootstrap resamples)}
\usage{
score_shd(boot.adj, alpha = 1, threshold=0, max.step = 500, 
blacklist = NULL, whitelist = NULL, print = FALSE)
}
\arguments{
  \item{boot.adj}{A p by p by B array, where B is the number of DAGs to be aggregated. 
  It records the adjacency matrices. It may be  the output of the "score" function.}
  \item{alpha}{a positive scalar: 
  alpha defines which member of the gSHD family should be used to aggregate the DAGs. 
  In general, the larger the alpha, the more aggressive of the aggregation, 
  in that less edges are retained leading to smaller FDR and less power, default  = 1}
  \item{threshold}{a scalar: it defines the frequency cut-off value, default =0, 
  meaning not to include any edge with (generalized) bootstrap selection frequency <0.5.}
    \item{max.step}{an integer: the maximum number of search steps, default = 500}
   \item{blacklist}{a p by p 0-1 matrix: if the (i,j)th-entry is "1", 
   then the edge i-->j will be excluded from the DAG during the search, default=NULL}
  \item{whitelist}{a p by p 0-1 matrix: if the (i,j)th-entry is "1", 
  then the edge i-->j will always be included in the DAG during the search, default=NULL}
 \item{print}{logic: whether print the step information, default= FALSE}
}
\details{
    score_shd aggregates an ensemble of DAGs to get a DAG 
    that minimizes the overall distance  to the ensemble. 
    A (generalized) structural hamming distance is 
    used as the distance metric on the space of DAGs of a given set of nodes. 
}
\value{
    A list with three components
    \item{adj.matrix}{adjacency matrix of the  aggregated DAG}
     \item{final.step}{a number recording how many search steps 
     are conducted before the procedure stops}
    \item{movement}{a matrix recording the selected operation, 
    addition, deletion or reversal of an edge, at each search step} 
}
\references{
Wang, R. and Peng, J. (2013) Learning directed acyclic graphs via bootstrap aggregating. arXiv:1406.2098
}
\author{
Wang, R and Peng, J.
}
\examples{
# read in data and true network
data(example)
Y.n=example$Y # data matrix 
true.dir=example$true.dir  #adjacency matrix of the data generating DAG
true.ske=example$true.ske  # skeleton graph of the data generating DAG
  

## (i) DAG learning  by hill climbing: no aggregation

# learn DAG using "BIC" score
temp=score(Y=Y.n, n.boot=0, score.type="BIC") 
adj=temp$adj.matrix

# Find the total number of skeleton edges 
# and the number of correct skeleton edges of the estimated DAG     
    tt=adj+t(adj) ##symmetrization
    correct.c=sum((tt>0)&(true.ske>0))/2    
    total.c=sum(tt>0)/2
    total.true=sum(true.ske>0)/2
    c(total.c, correct.c, total.true) ##326, 96, 109
    
## (ii) DAG learning by bagging   
set.seed(1)

# get the ensemble of DAGs from bootstrap resamples
temp.boot=score(Y.n, n.boot=10, score.type="BIC") 
boot.adj=temp.boot$adj.matrix

temp.bag=score_shd(boot.adj, alpha = 1) ##aggregating 
adj.bag=temp.bag$adj.matrix

# Find the total number of skeleton edges and the number of 
# correct skeleton edges of the estimated DAG     
 tt=adj.bag+t(adj.bag) ##symmetrization
  correct.c=sum((tt>0)&(true.ske>0))/2    
  total.c=sum(tt>0)/2
  total.true=sum(true.ske>0)/2
  c(total.c, correct.c, total.true) ##121, 89, 109 
 #Note: much less number of false edges and slightly less number of correct edges 
 #compared with the non-bagging result


## try different aggregation parameter alpha
#(i) alpha=0.5: less aggressive, more edges retained 
temp.bag=score_shd(boot.adj, alpha = 0.5) ##aggregating 
adj.bag05=temp.bag$adj.matrix

 # Find the total number of skeleton edges 
 # and the number of correct skeleton edges of the estimated DAG     
 tt=adj.bag05+t(adj.bag05) ##symmetrization
  correct.c=sum((tt>0)&(true.ske>0))/2    
  total.c=sum(tt>0)/2
  total.true=sum(true.ske>0)/2
  c(total.c, correct.c, total.true) ##132, 91, 109

#(ii) alpha=1.5: more aggressive, less edges retained 
temp.bag=score_shd(boot.adj, alpha = 1.5) ##aggregating 
adj.bag15=temp.bag$adj.matrix

 # Find the total number of skeleton edges 
 # and the number of correct skeleton edges of the estimated DAG     
 tt=adj.bag15+t(adj.bag15) ##symmetrization
  correct.c=sum((tt>0)&(true.ske>0))/2    
  total.c=sum(tt>0)/2
  total.true=sum(true.ske>0)/2
  c(total.c, correct.c, total.true) ##115, 87, 109
  
  #(iii) alpha=2: more aggressive, less edges retained 
  temp.bag=score_shd(boot.adj, alpha = 2) ##aggregating 
  adj.bag2=temp.bag$adj.matrix
  
   # Find the total number of skeleton edges 
   # and the number of correct skeleton edges of the estimated DAG     
   tt=adj.bag2+t(adj.bag2) ##symmetrization
    correct.c=sum((tt>0)&(true.ske>0))/2    
    total.c=sum(tt>0)/2
    total.true=sum(true.ske>0)/2
  c(total.c, correct.c, total.true) ##86, 67, 109

}