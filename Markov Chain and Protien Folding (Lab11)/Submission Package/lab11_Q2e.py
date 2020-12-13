# Author: Utkarsh Mali
# Adapted from Prof Nicolas

import numpy as np
import matplotlib.pyplot as plt
import random as random

# Combine function in to defenition.
print('[STATUS] Initializing...')

random.seed(10)


def Markov_Simluation(temp):
    # Initialize constants
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
    T = temp
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
            print(f"[STATUS] Markov Chain Model T= {temp}... {round(100 * step / steps, 1)}%")

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

    print("[STATUS] Simulated!\n")
    # Make simple plot
    plt.plot(magnet, label=f"T = {T}")

    # Return grid for imshow()
    return dipole_arr


# Plot figure
plt.figure()

# Save grids for imshow later on
# Plot change in magnetization over time
T1 = Markov_Simluation(1.0)
T2 = Markov_Simluation(2.0)
T3 = Markov_Simluation(3.0)
plt.xlabel("Time step")
plt.ylabel("Magnetization")
plt.title("Magnetization as a function of time")
plt.legend()
plt.savefig("Q2ept1.pdf")

# Plot final change in spin grid structure
f, ax = plt.subplots(1, 3)
f.set_figheight(3)
f.set_figwidth(7)

# Save data
data = [T1, T2, T3]

# Plot figures
for i in [0, 1, 2]:
    ax[i].imshow(data[i])
    ax[i].set_title(f'T = {i + 1.0}')
    ax[i].set_xlabel("Atoms (N)")
    if i == 0:
        ax[i].set_ylabel("Atoms (N)")

# Clean up plot
f.suptitle('Plots of Dipole Spin for Varying Temperatures', fontsize=16)
plt.savefig("Q2ept2.pdf")

print("[STATUS] All done!")
