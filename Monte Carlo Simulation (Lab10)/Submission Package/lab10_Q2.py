"""
author: Aslesha
"""

import numpy as np
from numpy.random import uniform


def f(r):
    """
    Return 1 if the squared norm of r is less than or equal to 1 and 0 
    otherwise. 
    """
    # generalize the check x^2+y^2 <= 1 to n dim
    squared_norm = np.dot(r, r) 
    if squared_norm <= 1:
        return 1
    
    # otherwise return 0
    return 0


# number of points
N = int(1e6)
# length of the hypercube
L = 2
# number of dimensions
D = 10
# volume of hypercube in D dim  
V = L**D

# to store sum of f and f^2
sum_f = 0
sum_fsq = 0

for i in range(N):
    # get a D dimensional array with points sampled uniformly at random 
    r = uniform(-1, 1, D)
    
    # accumulate the sum
    sum_f += f(r)
    sum_fsq += f(r)**2

# Get the final answer, generalization of eqn 10.33
I = V/N * sum_f

# calculate err
mean_f = 1/N * sum_f
mean_fsq = 1/N * sum_fsq
# find variance of f
var_f = mean_fsq - mean_f**2
# err given by eqn 10.32
err = V * np.sqrt(var_f) / np.sqrt(N)

print("The estimated volume of a 10-dimensional hypersphere is {:.2f} with the error of {:.2e}.".format(I, err))

