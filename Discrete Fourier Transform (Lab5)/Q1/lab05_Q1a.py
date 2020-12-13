# Author: Utkarsh Mali

import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
from tqdm import tqdm


def dft(y):
    N = len(y)
    c = np.zeros(N // 2 + 1, complex)
    for k in tqdm(range(N // 2 + 1)):
        for n in range(N):
            c[k] += y[n] * np.exp(-2j * np.pi * k * n / N)
    return c


# Loading data
data = np.loadtxt("sunspots.txt")

# Convert data into readable format
sunspots = data[:, 1]
month = data[:, 0]

# Question a
plt.figure()

# Make plot
plt.plot(month, sunspots, color='salmon')
plt.xlabel("Months since Jan 1749 ($t$)")
plt.ylabel("Number of Sunspots ($N$)")
plt.title("Number of sunspots plotted as a function of time")
plt.savefig("fig1.pdf")

# Question b
# Apply fft
sunspots_fft = dft(sunspots)
# square to get magnitude squared
ck_squared = abs(sunspots_fft) ** 2

# Make plot
plt.figure()
plt.plot(month[5:100], ck_squared[5:100])
plt.grid()
plt.xlabel("$k$")
plt.ylabel("Fourier coefficients $|c_k|^2$")
plt.title("Fourier coefficients $|c_k|^2$ as a function of $k$")
plt.savefig("fig2.pdf")

# Question c
# Find k after 0, add 10 since you clip by 10
k_index = np.argmax(ck_squared[10:]) + 10
print("The maximum value of |ck|^2 corresponds to k:", {k_index})

# Find period
period = len(sunspots[10:])/k_index
print(f"The period corresponding to the length of the cycle we found is: {round(period)}")
