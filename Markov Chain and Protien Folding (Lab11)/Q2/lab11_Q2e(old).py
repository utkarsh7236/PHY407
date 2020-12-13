# Author: Utkarsh Mali
# Adapted from Prof Nicolas

import numpy as np
import matplotlib.pyplot as plt
import random as random

print('[STATUS] Initializing...')
steps = int(1e5)
N = 20
J = 1.0

dipole_arr = np.empty((N, N), int)
for i in range(N):
    for j in range(N):
        dipole_arr[i, j] = random.choice([1, -1])


def energy_function(arr, J_=J):
    row, col = np.zeros([1, N], float), np.zeros(N, float)
    for i in range(0, N - 1):
        col += -J_ * arr[:, i] * arr[:, i + 1]
        row += -J_ * arr[i, :] * arr[i + 1, :]
    ret = sum(row[0]) + sum(col)
    return ret


def _acceptance_helper(beta, dE):
    return np.exp(- beta * dE)


def acceptance(beta, dE):
    if dE >= 0:
        P = _acceptance_helper(beta, dE)
    elif dE < 0:
        return True
    else:
        raise ValueError

    return False if random.random() > P else True


# Part e
T = 1.0
J = 1.0
kB = 1.0
beta = 1 / (kB * T)
energy = []
magnet = []

print("[STATUS] Running Monte-Carlo Simulation...")
for step in range(steps):
    if step % 25000 == 0:
        print(f"[STATUS] Markov Chain Model... {round(100 * step / steps, 1)}%")

    i, j = random.randint(0, N - 1), random.randint(0, N - 1)

    Eold = energy_function(dipole_arr)

    dipole_arr[i, j] *= -1

    Enew = energy_function(dipole_arr)

    dE = Enew - Eold

    if acceptance(beta, dE):
        Eold = Enew
        energy.append(Enew)

    else:
        dipole_arr[i, j] *= -1

    magnet.append(np.sum(dipole_arr))

print("\n[STATUS] Simulated!")
print("[STATUS] Plotting Figures...")

plt.figure()
plt.plot(magnet, label=f'T = {T}')

dipole_arr = np.empty((N, N), int)
for i in range(N):
    for j in range(N):
        dipole_arr[i, j] = random.choice([1, -1])

# Part d
T = 2.0
energy = []
magnet = []

print("[STATUS] Running Monte-Carlo Simulation...")
for step in range(steps):
    if step % 25000 == 0:
        print(f"[STATUS] Markov Chain Model... {round(100 * step / steps, 1)}%")

    i, j = random.randint(0, N - 1), random.randint(0, N - 1)

    Eold = energy_function(dipole_arr)

    dipole_arr[i, j] *= -1

    Enew = energy_function(dipole_arr)

    dE = Enew - Eold

    if acceptance(beta, dE):
        Eold = Enew
        energy.append(Enew)

    else:
        dipole_arr[i, j] *= -1

    magnet.append(np.sum(dipole_arr))

print("\n[STATUS] Simulated!")
print("[STATUS] Plotting Figures...")

plt.plot(magnet, label=f'T = {T}')

dipole_arr = np.empty((N, N), int)
for i in range(N):
    for j in range(N):
        dipole_arr[i, j] = random.choice([1, -1])

# Part e
T = 3.0
energy = []
magnet = []

print("[STATUS] Running Monte-Carlo Simulation...")
for step in range(steps):
    if step % 25000 == 0:
        print(f"[STATUS] Markov Chain Model... {round(100 * step / steps, 1)}%")

    i, j = random.randint(0, N - 1), random.randint(0, N - 1)

    Eold = energy_function(dipole_arr)

    dipole_arr[i, j] *= -1

    Enew = energy_function(dipole_arr)

    dE = Enew - Eold

    if acceptance(beta, dE):
        Eold = Enew
        energy.append(Enew)

    else:
        dipole_arr[i, j] *= -1

    magnet.append(np.sum(dipole_arr))

print("\n[STATUS] Simulated!")
print("[STATUS] Plotting Figures...")

plt.plot(magnet, label=f'T = {T}')

plt.show()
plt.legend()
print("[STATUS] All done!")
