import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft, irfft
from dcst import *

# Question a
data = np.loadtxt("dow2.txt")

# Create plot
plt.figure()
plt.plot(data, label="No Clipping")

# Apply fourier transform
fft_data = rfft(data)

# Clip to 2%
fft_data[int(0.02 * len(fft_data)):] = 0

# Apply reverse fourier transform
reverse_data = irfft(fft_data)

# Plot data
plt.plot(reverse_data, label="2% Clipping (DFT)")
plt.xlabel("Time ($t$)")
plt.ylabel("Index Valuation $\$$USD")
plt.title("DOW Index from 2004-2008")

# Question b
# Apply fourier transform
fft_data = dct(data)

# Clip to 2%
fft_data[int(0.02 * len(fft_data)):] = 0

# Apply reverse fourier transform
reverse_data = idct(fft_data)

plt.plot(reverse_data, label="2% Clipping (DCT)")
plt.legend()
plt.grid()

plt.savefig("fig4.pdf")
plt.show()
