# Author: Utkarsh Mali

import numpy as np
import matplotlib.pyplot as plt
from math import factorial, sqrt
from scipy.special import eval_hermite
from scipy.misc import derivative
from numpy import diff
from lab03_Q2_functions import gaussxw, gaussxwab

print("[STATUS] Initializing")


# Question 2a
# Define hermite polynomial function
def H(x, n):
    if n < 0:
        return 0
    if n == 0:
        return 1.
    elif n == 1:
        return 2. * x
    else:
        return 2 * x * H(x, n - 1) - 2 * (n - 1) * H(x, n - 2)


# set x values
x_arr = np.linspace(-4, 4, 100)


# define wavefunction
def wavefunction(n, x):
    c = 1 / (2 ** n * factorial(n) * np.sqrt(np.pi))
    return np.sqrt(float(c)) * np.exp((-x ** 2) / 2) * H(x, n)
    # return np.sqrt(float(c)) * np.exp((-x ** 2)/2) * eval_hermite(n,x)


print("[STATUS] Running A...")
# for loop over n
for n in [0, 1, 2, 3]:

    # f(x)
    fx = wavefunction(n, x_arr)

    # plot calculated f(x)
    plt.plot(x_arr, fx, label=f"n = {n}")

# labeling axis and making pretty
plt.minorticks_on()
plt.grid(color = 'grey',which = 'minor',linestyle = ":",linewidth = '0.1')
plt.grid(color = 'black',which = 'major')
plt.legend()
plt.xlim(min(x_arr), max(x_arr))
plt.xlabel("Position ($x$)")
plt.ylabel("Energy")
plt.title("Wavefunctions using hermite polynomials for n = [0,3]")
plt.savefig("fig1.pdf")

print("[STATUS] A complete running B...")
# Question 2b
plt.figure()
# new x
x_arr = np.linspace(-10, 10, 300)

# for loop over new n
n = 30

# for loop over x
fx = wavefunction(n, x_arr)

# plot calculated f(x)
plt.plot(x_arr, fx, label=f"n = {n}")

# labeling axis and making pretty
plt.minorticks_on()
plt.grid(color = 'grey',which = 'minor',linestyle = ":",linewidth = '0.1')
plt.grid(color = 'black',which = 'major')
plt.legend()
plt.xlim(min(x_arr), max(x_arr))
plt.xlabel("Position ($x$)")
plt.ylabel("Energy")
plt.title("Wavefunctions using hermite polynomials for n = 30")
plt.savefig("fig2.pdf")

print("[STATUS] B complete running C...")

# Question 2c
# initialize for loop

sums = []
for n in range(16):

    # def f(x) - position uncertainty with change for variable for infinity
    def f(x):
        z = x / (1 - x ** 2)
        integrand = abs(wavefunction(n, z)) ** 2
        return float(z ** 2 * integrand * (1 + x ** 2) / ((1 - x ** 2) ** 2))


    # set integration parameter and boundaries
    N = 100
    a = -1
    b = 1

    # apply gaussian quadrature
    x, w = gaussxw(N)
    s = 0.

    # 'integrate' over the sum
    for k in range(N):
        s += w[k] * f(x[k])

    sums.append(s)

# initialize sum
sums2 = []


# define partial derivative for wavefunction
def partial_derivative(n, x):
    c = float(1 / sqrt(2 ** n * factorial(n) * np.sqrt(np.pi)))
    e = np.exp((-x ** 2))
    hermite = (-x * H(x, n) + 2 * n * H(x, n - 1))
    return c * e * hermite


# next for loop
for n in range(16):

    # def f2(x) - momentum uncertainty with change for variable for infinity
    def f2(x):
        z = x / (1 - x ** 2)
        integrand = abs(partial_derivative(n, z)) ** 2
        return float(integrand * (1 + x ** 2) / ((1 - x ** 2) ** 2))


    # set integration parameter and boundaries
    N = 100
    a = -1
    b = 1

    # apply gaussian quadrature
    x, w = gaussxw(N)
    s = 0.

    # 'integrate' over the sum
    for k in range(N):
        s += w[k] * f2(x[k])

    sums2.append(s)

rms_position = list(np.sqrt(sums))
rms_momentum = list(np.sqrt(sums2))
energy = list(0.5 * (np.array(sums) + np.array(sums2)))

print('\nposition squared expectations (n=0,..,15):\n', np.array(sums), "\n")
print('momentum squared expectations (n=0,..,15):\n', np.array(sums2), "\n")
print('energy values calculated (n=0,..,15):\n', np.array(energy), "\n")
print("root mean squared position uncertainty:\n", np.array(rms_position), "\n")
print("root mean squared momentum uncertainty:\n", np.array(rms_momentum), "\n")
print("\n Position uncertainty at n = 5, ", round(rms_position[5], 5), "\n")
print('####################\n[STATUS] Finished\n####################')
