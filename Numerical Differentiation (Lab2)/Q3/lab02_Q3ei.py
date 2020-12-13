# Author: Utkarsh Mali
import numpy as np  # Package of choice
import math  # Required Package
import matplotlib.pyplot as plt
from tqdm import tqdm

# Question 1ei
# slit seperation
d_0 = 20. * 1.e-6

# Define transmission function q(u)
def q(u, d=d_0):
    alpha = np.pi / d
    beta = 0.5 * alpha
    return (np.sin(alpha * u) ** 2) * (np.sin(beta * u) ** 2)


# wavelength
lamda_0 = 500 * 1.e-9

# focal length
focal_0 = 1.

# total screen width
D_0 = 0.1

# numbers of slits
n = 10

# total width
w = n * d_0


# Defining I(x)
def I(x, u, focal=focal_0, lamda=lamda_0):  # du element of [-w/2, w/2]
    i = complex(0, 1)
    return np.sqrt(q(u)) * np.exp(i * 2 * np.pi * x * u / lamda * focal)


# Simpson's method of integration
# Note: I will be using the extended Simson's rule

# lower bound
a = -w / 2

# upper bound
b = w / 2

# Number of cuts
N = 1000

# step
h = w / N

x = np.linspace(-D_0 / 2, D_0 / 2, N)
Ix = []
for x_step in tqdm(x):
    # initialize even sum
    even = 0

    # initialize odd sum
    odd = 0
    # Odd terms
    for k in range(1, N, 2):
        odd += I(x_step, a + k * h)

    # Even terms
    for k in range(2, N, 2):
        even += I(x_step, a + k * h)

    # Complete final sum for integral
    Integral = (1. / 3) * h * (I(x_step, a) + I(x_step, b) + 4 * odd + 2 * even)
    Ix.append(abs(Integral) ** 2)

# Plotting result as fancy graph including heat map
fig, (ax, ax2) = plt.subplots(nrows=2)
plt.title("Intensity Heat Map and Graph with new $q(u)$")
x = np.array(x)
y = np.array(Ix)

# Creating heatmap
ex = [x[0] - (x[1] - x[0]) / 2., x[-1] + (x[1] - x[0]) / 2., 0, 1]
ax.imshow(y[np.newaxis, :], cmap="plasma", aspect="auto", extent=ex)
ax.set_xlim(ex[0], ex[1])
ax2.plot(x, y, color='salmon')

# Label axis
ax.set_ylabel("Intensity $I(x)$")
ax.set_xlabel("Position (m)")
ax2.set_ylabel("Intensity $I(x)$")
ax2.set_xlabel("Position (m)")
plt.tight_layout()
plt.savefig('fig4.pdf')
