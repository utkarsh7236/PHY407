# Author: Utkarsh Mali

import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft, irfft

# Question a
# Loading data
data = np.loadtxt("dow.txt")

# Plot data
plt.figure()
plt.plot(data, label="No Clipping")
plt.xlabel("Time ($t$)")
plt.ylabel("Index Valuation $\$$USD")
plt.title("DOW Index from 2006-2010")

# Question b
# Apply fourier transform
fft_data = rfft(data)

# Question c
# Clip to 10%
fft_data[int(0.1 * len(fft_data)):] = 0

# Question d
# Apply reverse fourier transform
reverse_data = irfft(fft_data)

# Plot data
plt.plot(reverse_data, "--", label="10% Clipping", linewidth=2)

# Question e
# Apply fourier transform
fft_data = rfft(data)

# Clip to 2%
fft_data[int(0.02 * len(fft_data)):] = 0

# Apply reverse fourier transform
reverse_data = irfft(fft_data)

# Plot data
plt.plot(reverse_data, label="2% Clipping")

# Savefig
plt.grid()
plt.legend()
plt.savefig("fig3.pdf")
