# LOGGLE (LOcal Group Graphical Lasso Estimation)

## Description
The R package `loggle` provides a set of methods that learn time-varying graphical models based on data measured over a temporal grid. `loggle` is motivated by the needs to describe and understand evolving interacting relationships among a set of random variables in many real applications, for instance, the gene regulatory networks over the course of organismal development, and the dynamic relationships between individuals in a community over a few years. `loggle` estimates time-varying graphical models under the assumption that the graph topology changes gradually over time.

`loggle` has been applied to S&P 500 stock price dataset, where the interacting relationships among stocks and among industrial sectors in a time period that covers the recent global financial crisis can be revealed. Detailed description of S&P 500 stock price dataset is in `?stockdata`.

`loggle` is available on [CRAN](https://CRAN.R-project.org/package=loggle) now. For more details on estimating time-varying graphical models and the package, please refer to: Yang, J. & Peng, J. (2018), **Estimating Time-Varying Graphical Models**, [arXiv:1804.03811](https://arxiv.org/abs/1804.03811). For codes and data used in the simulation and real data application in this paper, please refer to: [https://github.com/jlyang1990/loggle_test](https://github.com/jlyang1990/loggle_test).

## Dependencies
Please make sure to install the following package dependencies before using R package `loggle`. R with version later than 3.0.2 is needed.
```r
install.packages(c("Matrix", "doParallel", "igraph", "glasso", "sm"))
```

## Installation
The R package `loggle` can be installed from source files in the GitHub repository (R package `devtools` is needed):
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
