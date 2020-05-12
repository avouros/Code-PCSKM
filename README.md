This repository contains the code used to produce the results of the manuscript: [A semi-supervised sparse K-Means algorithm](https://arxiv.org/abs/2003.06973).

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
