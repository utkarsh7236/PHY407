import numpy as np
import matplotlib.pyplot as plt
from lab03_Q1_functions import trapezoidal_integral, simpson_integral
from lab03_Q1_gaussxw import gaussxwab
from tabulate import tabulate
from scipy.special import dawsn


########################## Reusing code from lab02 ##########################
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
                    Trapezoidal, Simpson's rule or Gaussian Quadrature.
    N : number of slices to use in the integration.
    x : Point to evaluate the integral at.
    """
    # use helper functions f(t) to evaluate dawson function
    return  np.e**(-x**(2)) * integral_func(N, f, 0, x)

#############################################################################


def gauss_wrapper(N, f, a, b):
    """
    Wrapper function for Gaussian Quardature to be compatible with 
    evaluate_dawson.
    """
    # get the weights and sample points
    x, w = gaussxwab(N, a, b)
    # the loop of summing wp[k] * f(xp[k]) implemented using numpy
    return np.sum(np.array(w) * f(np.array(x)))


def main():
    # evaluate at x = 4
    x = 4
    # store results as [[N], [Trapezoidal], [Simpson], [Gauss]]
    result = [[], [], [], []]
    for i in range(3, 12):
        N = 2 ** i
        # evaluate Dawson for various N using each integration method and store result
        result[0].append(N)
        result[1].append(evaluate_dawson(trapezoidal_integral, N, x))
        result[2].append(evaluate_dawson(simpson_integral, N, x))
        result[3].append(evaluate_dawson(gauss_wrapper, N, x))

    # title for the table        
    headers = ["N", "Trapezoidal rule", "Simpson's rule", "Gaussian Quadrature"]
    # create and print table using result
    tab = tabulate(np.array(result).T, headers=headers, tablefmt="pretty")
    print(tab)
    
    # the practical error estimation I_2N - I_N using np.diff as result already 
    # contains evaluation for I_N, I_2N, I_4N...
    error = np.diff(result)
    error[1] = 1 / 3 * error[0] # Need to multiply by 1/3 for Trapezoidal method
    error[2] = 1 / 15 * error[1] # Need to multiply by 1/15 for Simpson's method
    
    # plot the error estimation
    plt.figure()
    # since diff subtracts the conscutive entries of array, plot by excluding the 
    # first N as result[0][1:]
    plt.loglog(result[0][1:], error[1], label = "Trapezoidal")
    plt.loglog(result[0][1:], error[2], label="Simpson's")
    plt.loglog(result[0][1:], error[3], label = "Gaussian")
    plt.legend()
    plt.xlabel("N")
    plt.ylabel("Error")
    plt.title("Practical error estimation for different integration methods")
    plt.savefig("Q1ai.pdf")


    # evaluating Dawson function using scipy.special.dawsn
    dawsn_eval = dawsn(4)
    # plot the absolute relative error by taking dawsn as the True value
    plt.figure()
    # log-log plot for absolute relative error for all 3 methods
    plt.loglog(result[0], np.abs(dawsn_eval - result[1]) / dawsn_eval, label = "Trapezoidal")
    plt.loglog(result[0], np.abs(dawsn_eval - result[2]) / dawsn_eval, label = "Simpson's")
    plt.loglog(result[0], np.abs(dawsn_eval - result[3]) / dawsn_eval, label = "Gaussian")
    plt.legend()
    plt.xlabel("N")
    plt.ylabel("Error")
    plt.title("Absolute relative error using scipy.special.dawsn as the true value")
    plt.savefig("Q1aii.pdf")
    

if __name__ == "__main__":
    main()
