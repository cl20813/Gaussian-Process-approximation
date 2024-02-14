# Gaussian Process approximation.

The goal of this project is to approximate likelihood by using subset of data which follows gaussian process.

## Exercise 1: Simulate data using exponential kernel and try to fit variance only with matern kernel and do the opposite.
  1. If we try to fit matern when exponential kernel is true, fitted variance will go to infinity to compensate matern's faster decaying speed.
  2. On the other hand, if we try to fit exponential kernel when true model is mater, fitted variance will go to 0 to compensate exp's slower decaying speed.

  -[Experiment result](https://github.com/cl20813/Gaussian-Process-approximation/blob/567fc1abc0f2e12e7582635b54813c3ec11268d6/Exercises/Fit%20matern_true%20exp.pdf)


Click the link to review expectation and variance of theta(variance).
[Click here](https://stats.stackexchange.com/questions/427332/variance-of-quadratic-form-for-multivariate-normal-distribution)


  
## Exercise 2: Using eigen-values and eigen-vectors to assess the quality of the model.

  -[Experiment result](Exercises/Diagnostics_of_covariance_matrix_using_eigenvalue.ipynb)
  

## Exercise 3: Extension of exercise 2 to spatio-temporal covariance matrices in Stein(2004).

  -[Experiment result](Exercises/Spat_tmp_cov_exercise_stein_2004_python.ipynb)
  

## Exercise 4: Now we investigate the behavior of (eigen vector x complex exponential ) 

  -[Explanation](Exercises/Isometry_same_norm_spectrum_1_12.pdf)
  See the chapter 2.6 "Two corresponding Hilbert spaces" in Stein's book. Note that there is no randomness when working with spectral density.

  -[Python code](Exercises/Isometry_same_norm_spectrum_1_12.ipynb)

## Exercise 5: Now we investigate the behavior of (eigen vector x complex exponential ) 

1. Exercise 5-1: Observe the difference of isotrpoic model and anisotropic model.
2. Exercise 5-2: Let $f_1(w)= o(w^2)$ and $f_2(w)= o(w^4)$. If we look at $\biggl| \sum_j v_{jk} e^{i \omega x_j} \biggr|$, $f_2$ is lower than $f_1$ at w=0 and the difference decreases as w goes to high frequency. Here we use BLP for linear coefficients, $argmin_c  cov(  Z_0 - C^\top Z,  Z_0 - C^\top Z )   = \hat{C}= \Gamma^{-1} \gamma.$

3. Exercise 5-3: Observe that most of the variance of BLP is at high frequencies, another reason why we should focus on local behaviors.
   See "Predicting random fields", Stein (1999).

   -[Python code](Exercises/Experiment5.ipynb)

## Exercise 6. Correlation and relative standard errors of empirical semi-variogram. 

  -[Explanation](Exercises/Experiment6_acf_semivariogram.pdf)
  
  -[Python code](Exercises/Experiment6_acf_semivarogram.ipynb)
  
  
  


 


