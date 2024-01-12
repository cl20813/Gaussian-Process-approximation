# Gaussian Process approximation.

The goal of this project is to approximate likelihood by using subset of data which follows gaussian process.

## Exercise 1: Simulate data using exponential kernel and try to fit variance only with matern kernel and do the opposite.
  1. If we try to fit matern when exponential kernel is true, fitted variance will go to infinity to compensate matern's faster decaying speed.
  2. On the other hand, if we try to fit exponential kernel when true model is mater, fitted variance will go to 0 to compensate exp's slower decaying speed.

  -[Experiment result](https://github.com/cl20813/Gaussian-Process-approximation/blob/567fc1abc0f2e12e7582635b54813c3ec11268d6/Exercises/Fit%20matern_true%20exp.pdf)


Click the link to review expectation and variance of theta(variance).
[Click here](https://stats.stackexchange.com/questions/427332/variance-of-quadratic-form-for-multivariate-normal-distribution)


  
## Exercise 2: Using eigen-values and eigen-vectors to assess the quality of the model.

  -[Experiment result](https://github.com/cl20813/Gaussian-Process-approximation/blob/96833cdb8ed31c74473dfea7213204bed6111942/Exercises/Diagnostics%20of%20models%20using%20eigenvalues.pdf)
  

## Exercise 3: Extension of exercise 2 to spatio-temporal covariance matrices in Stein(2004).

  -[Experiment result](https://github.com/cl20813/Gaussian-Process-approximation/blob/main/Exercises/Spat_tmp_cov_exercise_stein_2004_python.ipynb)
  

## Exercise 4: Extension of exercise 2 to spatio-temporal covariance matrices in Stein(2004).

  -[Experiment result](Exercises/Isometry_same_norm_spectrum_1_12.pdf)
  


 


