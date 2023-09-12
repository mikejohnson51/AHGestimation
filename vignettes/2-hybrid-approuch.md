# A Hybrid Method To Model Hydraulic Geometry Relations

Traditionally AHG relations were derived using Ordinary Least Squares (OLS) on the logarithmic transformation of the variables. [Studeis] (https://edx.hydrolearn.org/assets/courseware/v1/732bdf6048116238221bf6d282b23434/asset-v1:Utah_State_University+CEE6400+2019_Fall+type@asset+block/2020_09_08.pdf) have highlighted the potential problem of using this approach given that when the log transformed values are back-transformed, the estimates are effectively medians, and not means of the estimates, resulting in a general low bias. 

Nonlinear Least Square (NLS) regression provides more flexible curve fitting through an iterative optimization. NLS approaches require a specified starting value for each un-known parameter to ensure the solver converges on a global rather than a local minimum. When suboptimal starting values are provided, NLS solvers may converge on a local min-imum, or, not at all.

Both OLS and NLS minimizing the error between the predicted and observed values for each relation. Such fitting procedure does not grantee system achieving the continuity and therefore an objective function that considers flow continuty must be imposed. This software uses Evolutionary computing approaches (EA) provide stochastic search algorithms in-spired by the principles of natural selection in which the “fittest” individuals from each population pass on their parameters to the next population, mimicking processes like se-lection, crossover, and mutation. The notion of “fittest” can be prescribed based on a user-defined objective function (OF).

$OF = \sum[TW_{nRMSE}  ,V_{nRMSE}  ,Y_{nRMSE} ]$

The conventional approuches in literature would 

Therefore, when a time series for TW, V and Y are available, an objective function can be de-fined as the sum of each nRMSE and solved with an elitist version of a Ge-netic Algorithm (GA) using the [R package GA](https://www.jstatsoft.org/article/view/v053i04).

In the GA solver, solutions can be penalized for violating the continuity constraints by adding a weight to the resulting objective function. Using this type of objective function assumes the three nRMSE metrics are not competing, and that lowering, for example, TW nRMSE, will not increase nRMSE in the U or Y relation. To test whether fitting AHG is a multi-criterion optimization problem, the Non-dominated Sorting Genetic Algorithm II (NSGA-II) elitist multi-objective Genetic Algorithm was implemented to minimize the Multi-Criterion OF (MCOF) defined in Equation 12, subject to the constraints shown in Equations 12,13.

$MCOF = TW_{nRMSE}  ,V_{nRMSE}  ,Y_{nRMSE}$

$1-allowance ≤ b+f+ m$

$1-allowance ≤  a × c× k$
``` r
library(FHGestimation)
usgs_obs <- load("extra/usgs_obs_hydraulics.rda")

index = 50
  
tmp = usgs_obs[[index]] %>% 
    filter(as.Date(date) > as.Date('2010-01-01')) %>% 
  select(date,Q, TW, V, Y = Ymean)

calc_ga(tmp$Q, tmp$Y, tmp$V, 
            tmp$TW, allowance = .05, r, type = "nrmse") 
```

## The Automated Workflow

When three hydraulic states are provided, continuity is computed for the OLS and NLS solutions. If continuity for NLS or OLS is within a prescribed allowance (default = 0.05), the summary results are returned, and the solver concludes. If continuity is not met, the NSGA-II algorithm is run with 500 generations, a population size of 100, a crossover probability of 0.7, a crossover distribution index of 5, a mutation probability of 0.2, and a mutation distribution index of 10 using the mco R package 

![Flowchart](/man/figures/diagram.png)

``` r
(x = fhg_estimate(Q = tmp$Q, 
             V = tmp$V, 
             TW = tmp$TW, 
             Y = tmp$Y, 
             tmp = .05)$output)
```
