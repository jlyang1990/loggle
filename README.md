# LOGGLE (LOcal Group Graphical Lasso Estimation)

## Description
The R package `loggle` provides a set of methods that learn time-varying graphical models, where these models can be used to identify the dynamic relationships between a bunch of variables over a time period. `loggle` is motivated by the cases where the observations of variables are measured over a time grid, such that the direct interactions between these variables can evolve over time. For example, the gene regulatory network over the course of organismal development, and the dynamic relationship between individuals in a community over several years. `loggle` estimates the time-varying graphical models under smoothness assumptions on both the covariance matrix and graph structure. 

`loggle` has been applied to S&P 500 stock price dataset from 2007 to 2016, where the evolution pattern of relationships between stocks during that period can be revealed. Detailed description and analysis of S&P 500 stock price dataset is in `?stockdata`.

## Dependencies
Please make sure to install the following package dependencies before using R package `loggle`. R with version later than 3.0.2 is prerequisite.
```r
install.packages(c("Matrix", "doParallel", "igraph", "glasso", "sm"))
```

## Installation
The R package `loggle` can be installed from source files in the GitHub repository (installing R package `devtools` is prerequisite):
```r
library(devtools)
install_github(repo="jlyang1990/loggle")
```

## Main functions
* `loggle`: learn time-varying graphical models for a given set of tuning parameters.
* `loggle.cv`: conduct model selection via cross validation for learning time-varying graphical models.
* `loggle.cv.select`: conduct model selection for time-varying graphical models based on cross validation results from `loggle.cv`.
* `loggle.cv.vote`:  learn time-varying graphical models for a given set of tuning parameters via cv.vote.
* `loggle.refit`: conduct model refitting given learned time-varying graph structures.

## Contact
Please report any bugs to jlyang@ucdavis.edu.
