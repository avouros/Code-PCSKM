This repository contains the code used to produce the results of the manuscript: [A semi-supervised sparse K-Means algorithm](https://www.sciencedirect.com/science/article/abs/pii/S0167865520304268) ([arxiv version](https://arxiv.org/abs/2003.06973)).

# Code-PCSKM

__exeSimus.m:__ Runs the whole analysis and stores the results inside the _./GenRes/results_ folder. This file contains the following options:
 - __DETERM:__ 0/1 start without or with a random seed.
 
 - __JMPCKM_OVERLOAD:__ 0/1 use overloaded or non-overloaded MPCK-Means. The [WekaUT](http://www.cs.utexas.edu/users/ml/risc/code/) library is used for the MPCK-Means algorithm. See [Bilenko, M., et al. (2004)](https://dl.acm.org/doi/abs/10.1145/1015330.1015360).

- __CONSTR_PERC:__ 0/1 use a flat number of constraints or percentages based on size.

- __LOG:__ (0) no log file and no display, (1) log file only, (2) display only, (else) both display and log file.

- __constraints_type:__ Type of constraints to use; 0/1 to activate ML and/or CL, when both 1 then equal number of constriants per type is selected when either -1 then random constraints are picked from all the available constraints.

- __constraints_number:__ flat or percentage of constraints to use.

- __citer:__ number of iterations per constraints

- __sstep:__ sparsity parameter values to be tested form 1.1 to sqrt(dimensions) with step _sstep_.

- __maxIter:__ iterations for algorithm to reach convergence.

- __kfolds:__ selection of _k_ for k-fold validation.

__CVstatsPer.m:__ Generates statistics about the data sets such as percentage of used constraints during the k-fold validation.



## Citations for software and code that we have used in this project

**Density K-Means++:**

[Nidheesh, N., KA Abdul Nazeer, and P. M. Ameer. "An enhanced deterministic K-Means clustering algorithm for cancer subtype prediction from gene expression data." Computers in biology and medicine 91 (2017): 213-221.](https://www.sciencedirect.com/science/article/pii/S0010482517303402)

MATLAB code was based on the R implementation of the algorithm; code: [`dkmpp_0.1.0`](https://github.com/nidheesh-n/dkmpp)

**MPCK-Means:**

[Bilenko, Mikhail, Sugato Basu, and Raymond J. Mooney. "Integrating constraints and metric learning in semi-supervised clustering." Proceedings of the twenty-first international conference on Machine learning. 2004.](https://dl.acm.org/doi/abs/10.1145/1015330.1015360)

Modified [WekaUT](http://www.cs.utexas.edu/users/ml/risc/code/) in order to read initial centroids from text files and write results to text files.

**Sparse clustering:**

[Witten, Daniela M., and Robert Tibshirani. "A framework for feature selection in clustering." Journal of the American Statistical Association 105.490 (2010): 713-726.](https://amstat.tandfonline.com/doi/abs/10.1198/jasa.2010.tm09415)

[Brodinová, Šárka, et al. "Robust and sparse k-means clustering for high-dimensional data." Advances in Data Analysis and Classification (2017): 1-28.](https://link.springer.com/article/10.1007/s11634-019-00356-9)

MATLAB code was based on the R implementation of the algorithm; packages: [`sparcl`](https://cran.r-project.org/web/packages/sparcl/index.html) and [`wrsk`](https://github.com/brodsa/wrsk)


