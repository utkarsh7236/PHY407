# Author: Utkarsh Mali

import numpy as np
import matplotlib.pyplot as plt

# import constants
from scipy.constants import hbar, eV, angstrom, electron_mass, elementary_charge

# define constants
a_0 = 10  # in eV0)
width = 5 * angstrom


# Question 2b
# Define H(m,n)
def Hmatrix(m, n, L=width, a=a_0):
    # m == n case
    if m == n:
        # divide by eV for unit conversion
        num = (np.pi ** 2) * (hbar ** 2) / eV * (m ** 2)
        dem = 2 * electron_mass * (L ** 2)
        return 0.5 * a + num / dem

    # Both odd/even case
    elif (m != n) & (m % 2) == (n % 2):
        return 0

    # Either odd/even case
    elif m != n:
        num = 8 * a * m * n
        dem = (np.pi ** 2) * ((m ** 2) - (n ** 2)) ** 2
        return -num / dem


# Question 2c
# Create empty matrix
m, n = 10, 10
mmax = m
nmax = n
H = np.empty([m, n])

# compute the eigenvalues now
for m in range(1, mmax + 1):
    for n in range(1, nmax + 1):
        H[m - 1, n - 1] = Hmatrix(m, n)

# we use eigh since it is a symmetrical real matrix
E_val, E_vec = np.linalg.eigh(H)

# Display results
print("First Ground State Energy is:", round(E_val[0], 4), "eV")

# Question 2d
# Create empty matrix
m, n = 100, 100
mmax = m
nmax = n
H = np.empty([m, n])

# compute the eigenvalues now
for m in range(1, mmax + 1):
    for n in range(1, nmax + 1):
        H[m - 1, n - 1] = Hmatrix(m, n)

# we use eigh since it is a symmetrical real matrix
E_val_, E_vec_ = np.linalg.eigh(H)

# Display results
print("Comparison of eigenvalues for first 10 points:\n", E_val, "\n", E_val_[:10])


# Question 2e
# Define V(x)
def V(x_0, a=a_0, L=width):
    return (a * x_0) / L


# Define Psi(x)
def Psi(x, excitation, L=width):
    s = 0
    # sum over n for a specific psi
    for n in range(1, nmax):
        s += E_vec_[excitation + 1, n] * np.sin(n * np.pi * x / L)
    return s


# Set values for x to compute
x = np.linspace(0, width, 250)

# Empty arrays for Psi(x)
Psi1_arr = []
Psi2_arr = []
Psi3_arr = []


# Define probability density function
def probability_density(x, excitation):
    return abs(Psi(x, excitation)) ** 2


# loop over x for first few energies of psi
for x_step in x:
    Psi1_arr.append(probability_density(x, 0))
    Psi2_arr.append(probability_density(x, 1))
    Psi3_arr.append(probability_density(x, 2))

# define gaussian quadrature for integrating
from lab04_Q2_functions import gaussxw


# Find normalization constant by integrating using gaussian quadrature
def Normalize(excitation):
    # steps
    N = 100

    # start, stop
    a = 0
    b = width

    # apply gaussxw
    x, w = gaussxw(N)
    xp = 0.5 * (b - a) * x + 0.5 * (b + a)
    wp = 0.5 * (b - a) * w

    # initialize sum and loop
    A = 0
    for k in range(N):
        A += wp[k] * probability_density(xp[k], excitation)

    # return square root
    return float(np.sqrt(A))


# Normalize desired energy states
Normed_Psi1 = Normalize(0) * np.array(Psi1_arr[0])
Normed_Psi2 = Normalize(1) * np.array(Psi2_arr[0])
Normed_Psi3 = Normalize(2) * np.array(Psi3_arr[0])

# plot figure
plt.figure()

# add all plots onto same figure
plt.plot(x, Normed_Psi1, label="Ground State")
plt.plot(x, Normed_Psi2, label="First Excited State")
plt.plot(x, Normed_Psi3, label="Second Excited State")

# label axis
plt.xlabel("Position x ($m$)")
plt.ylabel("Normalized $|\psi (x)|^2$")

# title and legend
plt.title("Psi Squared as a function of x")
plt.legend()

# save/show
plt.savefig("fig1.pdf")
