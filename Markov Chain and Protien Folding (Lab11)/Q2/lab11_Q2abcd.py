# Author: Utkarsh Mali
# Adapted from Prof Nicolas

import numpy as np
import matplotlib.pyplot as plt
import random as random

# Initialize constants
print('[STATUS] Initializing...')
steps = int(1e5)
N = 20
J = 1.0

# Initialize dipole array to store spins
dipole_arr = np.empty((N, N), int)
for i in range(N):
    for j in range(N):
        dipole_arr[i, j] = random.choice([1, -1])


# Define energy function as required in part a
def energy_function(arr, J_=J):
    row, col = np.zeros([1, N], float), np.zeros(N, float)
    for i in range(0, N - 1):
        col += -J_ * arr[:, i] * arr[:, i + 1]
        row += -J_ * arr[i, :] * arr[i + 1, :]
    ret = sum(row[0]) + sum(col)
    return ret


# Define helper function for the acceptance
def _acceptance_helper(beta, dE):
    return np.exp(- beta * dE)


# Define acceptance function the same was as Prof Nicolas
def acceptance(beta, dE):
    if dE >= 0:
        P = _acceptance_helper(beta, dE)
    elif dE < 0:
        return True
    else:
        raise ValueError

    return False if random.random() > P else True


# Initialize constants for part b
T = 1.0
J = 1.0
kB = 1.0
beta = 1 / (kB * T)

# Initialize empty arrays
energy = []
magnet = []

# Start Markov Chain Model for simluation
print("[STATUS] Running Monte-Carlo Simulation...")
for step in range(steps):
    # Update counter
    if step % 25000 == 0:
        print(f"[STATUS] Markov Chain Model... {round(100 * step / steps, 1)}%")

    # Choose i and j at random
    i, j = random.randint(0, N - 1), random.randint(0, N - 1)

    # Compute old energy
    Eold = energy_function(dipole_arr)

    # Flip spin randomly
    dipole_arr[i, j] *= -1

    # Calculate energy with flipped spin
    Enew = energy_function(dipole_arr)

    # Calculate the change in energy
    dE = Enew - Eold

    # Run acceptance limitation
    if acceptance(beta, dE):
        # Set new E if accepted and update energy
        Eold = Enew
        energy.append(Enew)

    else:
        # Re-flip dipole if no accepted
        dipole_arr[i, j] *= -1

    # Finally update magnetic list.
    magnet.append(np.sum(dipole_arr))

print("\n[STATUS] Simulated!")
print("[STATUS] Plotting Figures...")

# Plot figure
plt.figure()
plt.plot(magnet)
plt.xlabel("Time step")
plt.ylabel("Magnetization")
plt.title("Magnetization as a function of time")
plt.savefig("Q2c.pdf")

print("[STATUS] All done!")
