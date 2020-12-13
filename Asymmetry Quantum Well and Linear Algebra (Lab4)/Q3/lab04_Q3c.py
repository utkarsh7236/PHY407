import numpy as np
import matplotlib.pyplot as plt

# Question 3c, Section a
# define error
error = 1e-6


def newton_method(f, derv_f, x1, epsilon):
    xn = x1
    while abs(f(xn)) > epsilon:
        fxn = f(xn)
        dfxn = derv_f(xn)
        xn = xn - fxn / dfxn
    return xn


# helper for calculating the midpoint
def _midpoint(x1, x2):
    return 0.5 * (x1 + x2)


# helper for calculating the same sign
def _same_sign(f_x1, f_midpoint):
    return np.sign(f_x1) == np.sign(f_midpoint)


# Defining the Binary Search Method
def binary_search(function, x1, x2, epsilon):
    # set initial delta
    delta = abs(x1 - x2)

    # keep looping until we have a certain accuracy
    while delta > epsilon:

        # compute midpoint and its function values at x1 and midpoint
        x_prime = _midpoint(x1, x2)
        f_x1 = function(x1)
        f_midpoint = function(x_prime)

        # check for same sign and change x1, x2 accordingly
        if _same_sign(f_x1, f_midpoint):
            x1 = x_prime
        else:
            x2 = x_prime

        # update delta
        delta = abs(x1 - x2)
    # return the next iteration
    return _midpoint(x1, x2)


# Defining function for binary search method
def func(x):
    return (5 * np.exp(-x)) + x - 5


# define second fucntion for relaxation method
def func2(x):
    return -(5 * np.exp(-x)) + 5


def derv_func(x):
    return (-5 * np.exp(-x)) + 1


# Initializing Newton's Methods
# set the initial value of x
x0 = 1

# initialize difference in x to 1
dx = 1

# list to store subsequent x0
relaxation = [x0]

while dx > error:
    # append the new value to the list
    relaxation.append(func2(relaxation[-1]))
    # set the dx to the difference of last two x0s
    dx = np.abs(relaxation[-1] - relaxation[-2])

print("Relaxation Method:", relaxation[-1])
print("Newton's Method", newton_method(func, derv_func, 4, error))
print("Binary Search Method:", binary_search(func, x1=1, x2=100, epsilon=error))

# Question3c, Section b
# import useful constants
from scipy.constants import h, c, Boltzmann

# define lambda
lamda = 502 * 10 ** (-9)

# Compute T
num = h * c
dem = Boltzmann * lamda * binary_search(func, x1=1, x2=100, epsilon=error)
T = num / dem
print("Surface Temperature of the Sun:", T)
