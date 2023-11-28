# Gaussian Process approximation.

The goal of this project is to approximate likelihood by using subset of data which follows gaussian process.

## Exercise 1: Simulate data using exponential kernel and try to fit variance only with matern kernel and do the opposite.
  1. If we try to fit matern when exponential kernel is true, fitted variance will go to infinity to compensate matern's faster decaying speed.
  2. On the other hand, if we try to fit exponential kernel when true model is mater, fitted variance will go to 0 to compensate exp's slower decaying speed.
     -  [Experiment result](Exercises/Diagnostics of models using eigenvalues.pdf)
     -  [R code](Exercises/Diagnostics of models using eigenvalues.Rmd)


 




See below for variance of theta
https://stats.stackexchange.com/questions/427332/variance-of-quadratic-form-for-multivariate-normal-distribution
