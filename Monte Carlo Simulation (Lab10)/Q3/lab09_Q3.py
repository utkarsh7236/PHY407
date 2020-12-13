"""
author: Aslesha
"""

import numpy as np
from numpy.random import uniform, rand, normal
import matplotlib.pyplot as plt
import matplotlib

# set the fontsize to 14 for the plots
font = {'size' : 14}
matplotlib.rc('font', **font)

############ Q3a ############
def fa(x):
    """
    Return the value of integrand eqn 4 at x.
    """
    return (x**(-1/2)) / (1 + np.e**x)

def pa(x):
    """
    Distribution to sample the points from.
    """
    return 1 / (2 * x**(1/2))

def transform_a(z):
    """
    Return the transformation of z in terms of p(x).
    z: random sample in [1,0)
    """
    return z**2


############ Q3b ############
def fb(x):
    """
    Return the value of integrand eqn 6 at x.
    """
    return np.exp(-2 * np.abs(x-5))

def pb(x):
    """
    Distribution to sample the points from.
    """
    return 1 / np.sqrt(2*np.pi) * np.exp(- (x-5)**2 / 2) 


def mean_val_method(a, b, f, N, hist_range, savename, part_a=True):
    """
    Compute the integral of f using mean value method and plot the histogram 
    of the estimation.
    a, b : ends 
    f : integrand function
    N : number of sample points
    hist_range: range to plot the histogram
    part_a : boolean for title purpose
    """
    # result of mean value method in 100 runs
    mean_val_result = []
    
    for j in range(100):
        # to store sum of f 
        sum_f = 0
        
        for i in range(N):
            # get points sampled uniformly at random 
            r = uniform(a, b)
            
            # accumulate the sum
            sum_f += f(r)
        
        # Get the final answer, eqn 10.30
        I = (b-a)/N * sum_f 
        
        # collect all the results
        mean_val_result.append(I)
    
    # plot the histogram
    plt.figure(figsize=(8,6))
    plt.hist(mean_val_result, 10, range=hist_range)
    plt.ylabel("Frequency")
    plt.xlabel("Estimated integration")
    if part_a:
        plt.title("Mean value method estimation of $\int_0^1 x^{-1/2}/(1 + e^x) dx$")
    else:
        plt.title("Mean value method estimation of $\int_0^{10} exp(-2|x-5|) dx$")
    plt.tight_layout()
    plt.savefig(savename)
    plt.show()
    

def imp_sampling_method(f, p, N, hist_range, savename, part_a=True):
    """
    Compute the integral of f using importance sampling method and plot the 
    histogram of the estimation.
    f : integrand function
    p : probability distribution function used for drawing the points 
    N : number of sample points
    hist_range: range to plot the histogram
    part_a : boolean for title purpose
    """
    # result of importance sampling method in 100 runs
    imp_samp_result = []
    
    for j in range(100):
        # to store sum of f 
        sum_f = 0
        
        for i in range(N):
            if part_a:
                # get a point at random 
                z = rand()
                # transform it to p(x) distribution
                x = transform_a(z)
            else:
                x = normal(5, 1)
            
            # accumulate the sum according to eqn 10.42 and 10.39
            sum_f += f(x) / p(x)
        
        # Get the final answer, eqn 10.42 and 10.39
        I = 1/N * sum_f 
        
        imp_samp_result.append(I)
    
    # plot the histogram
    plt.figure(figsize=(8,6))
    plt.hist(imp_samp_result, 10, range=hist_range)
    if part_a:
        plt.title("Importance sampling estimation of $\int_0^1 x^{-1/2}/(1 + e^x) dx$")
    else:
        plt.title("Importance sampling estimation of $\int_0^{10} exp(-2|x-5|) dx$")
    plt.ylabel("Frequency")
    plt.xlabel("Estimated integration")
    plt.tight_layout()
    plt.savefig(savename)
    plt.show()
    

def q3a():
    """
    Run q3a.
    """
    # boundary of integration
    a = 0
    b = 1
    # number of points for integration
    N = 10000    
    
    # run mean value method
    mean_val_method(a, b, fa, N, [0.8, 0.88], "q3a_mean.pdf")
    
    # run the importance sampling method
    imp_sampling_method(fa, pa, N, [0.8, 0.88], "q3a_samp.pdf")
    
    
def q3b():
    """
    Run q3b.
    """
    # boundary of integration
    a = 0
    b = 10
    # number of points for integration
    N = 10000    
    
    # run mean value method
    mean_val_method(a, b, fb, N, [0.93, 1.07], "q3b_mean.pdf", part_a=False)
    
    # run the importance sampling method
    imp_sampling_method(fb, pb, N, [0.93, 1.07], "q3b_samp.pdf", part_a=False)
    
    
# run all the experiments
q3a()
q3b()
    
    


