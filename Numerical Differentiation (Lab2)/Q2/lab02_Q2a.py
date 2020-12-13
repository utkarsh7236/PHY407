import numpy as np
from lab02_Q2_functions import trapezoidal_integral, simpson_integral
from scipy.special import dawsn
from time import perf_counter


def f(t):
    """
    Function inside the integral in the Dawson function i.e. e^(t^2).
    t : point to evaluate f at.
    """
    return np.e**(t**2)


def evaluate_dawson(integral_func, N, x):
    """
    Evaluate the Dawson function at x using integral_func integration for N
    slices.
    integral_func : Integration function to use for integral evaluation. 
                    Either Trapezoidal or Simpson's rule.
    N : number of slices to use in the integration.
    x : Point to evaluate the integral at.
    """
    # use helper functions f(t) to evaluate dawson function
    return  np.e**(-x**(2)) * integral_func(N, f, 0, x)
    

# Q2 (a) (i) Comparing the trapezoidal and Simpson's rule to scipy.special.dawsn

# Point to evaluate Dawson function at
x = 4
# Number of slices
N = 8

# evaluating Dawson function using trapezoidal integral 
trapezoidal_eval = evaluate_dawson(trapezoidal_integral, N, x)

# evaluating Dawson function using Simpson's integral 
simpson_eval = evaluate_dawson(simpson_integral, N, x)

# evaluating Dawson function using scipy.special.dawsn
dawsn_eval = dawsn(4)

# Print the results
print("The value of Dawson function at x={}".format(x))
print("\t using Trapezoidal integral is {}".format(trapezoidal_eval))
print("\t using Simpson's integral is {}".format(simpson_eval))
print("\t using scipy.special.dawsn is {} \n".format(dawsn_eval))


# Q2 (a) (ii) Finding the nuber of slices to get an error of O(10^{-9}) using the dumb way

# for Trapezoidal rule start from N = 2^3 = 8 as above
n = 3
trapezoidal_eval = evaluate_dawson(trapezoidal_integral, 2**n, x)
# use np.allclose with absolute tolerance of 10^{-9}
while not np.allclose(trapezoidal_eval, dawsn_eval, atol=1e-09):
    n += 1
    trapezoidal_eval = evaluate_dawson(trapezoidal_integral, 2**n, x)

print("For the Trapezoidal rule it takes 2^{} slices to approximate the integral"
      " with an absolute error of O(10^(-9)).".format(n))


# for Simpson's rule start from N = 2^3 = 8 as above
n = 3
simpson_eval = evaluate_dawson(simpson_integral, 2**n, x)
# use np.allclose with absolute tolerance of 10^{-9}
while not np.allclose(simpson_eval, dawsn_eval, atol=1e-09):
    n += 1
    simpson_eval = evaluate_dawson(simpson_integral, 2**n, x)
    

print("For the Simpson's rule it takes 2^{} slices to approximate the integral"
      " with an absolute error of O(10^(-9)).\n".format(n))


# Q2 (a) (ii) Timing the two methods using N=2^n slices found above

# For trapezoidal rule
# accumulate the time taken
total_time_trapezoidal = []
# n found from above to achieve an error of O(10^{-9}) 
n = 12
# number of times to run
num_iteration = 100
# run num_iteration times
for i in range(num_iteration):
    s = perf_counter()
    evaluate_dawson(trapezoidal_integral, 2**n, x)
    total_time_trapezoidal.append((perf_counter() - s) * 1000) # accumulate time in ms

print("For the Trapezoidal rule it takes an average of {:.3f}ms to approximate an"
      " integral with error of O(10^(-9)).".format(np.mean(total_time_trapezoidal)))


# For Simpson's rule
# accumulate the time taken
total_time_simpson = []
# n found from above to achieve an error of O(10^{-9}) 
n = 8
# number of times to run
num_iteration = 100
# run num_iteration times
for i in range(num_iteration):
    s = perf_counter()
    evaluate_dawson(simpson_integral, 2**n, x)
    total_time_simpson.append((perf_counter() - s) * 1000) # accumulate time in ms

print("For the Simpson's rule it takes an average of {:.3f}ms to approximate an"
      " integral with error of O(10^(-9)).".format(np.mean(total_time_simpson)))


# Q2 (a) (iii) Practical error estimation
# Use n1 = 32 and n2 = 64
N1 = 32
N2 = 64
# evaluate at x = 4
x=4

# Estimate dawson function using Trapezoidal rule for N1 and N2 steps at x = 4 
trapezoidal_I1 = evaluate_dawson(trapezoidal_integral, N1, x)
trapezoidal_I2 = evaluate_dawson(trapezoidal_integral, N2, x)
# Trapezoidal error estimation given by 1/3(I2 - I1)
trapezoidal_err_estimation = 1/3 * (trapezoidal_I2 - trapezoidal_I1)
print("The practical error estimation for the trapezoidal rule is {}".format(trapezoidal_err_estimation))


# Estimate dawson function using Simpson's rule for N1 and N2 steps at x = 4 
simpson_I1 = evaluate_dawson(simpson_integral, N1, x)
simpson_I2 = evaluate_dawson(simpson_integral, N2, x)
# Trapezoidal error estimation given by 1/15(I2 - I1)
simpson_err_estimation = 1/15 * (simpson_I2 - simpson_I1)
print("The practical error estimation for the simpson's rule is {}".format(simpson_err_estimation))



