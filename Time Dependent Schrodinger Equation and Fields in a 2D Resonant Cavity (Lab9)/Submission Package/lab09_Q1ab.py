# Author: Utkarsh Mali
# Solving the time-dependent Schrödinger equation with the Crank-Nicolson scheme

import matplotlib.pyplot as plt
import numpy as np

print("[STATUS] Initializing...")

# define constants as required in Q1a
h = 6.62607004e-34
hbar = h / (2 * np.pi)
L = 1e-8
m = 9.109e-31
sigma = L / 25
P = 1024
x_0 = L / 5
N = 3000
tau = 1e-18
a = L / P
x = np.linspace(-L / 2, L / 2, P)
p = np.arange(1, P, 1)
kappa = 500 / L


# Define a checker for normalization
def isPsiNormalized(psi):
    normed = np.sum(abs(psi ** 2)) * a
    return np.isclose(normed, 1)


# Define value of normalization constant
def normalizedPsiValue(psi):
    normed = np.sum(abs(psi) ** 2) * a
    return normed


# Define position expectation computer
def positionExpectation(psi, x):
    ret = np.sum(np.conj(psi) * x * psi) * a
    return ret


# Compute psi without psi0
xdiff = (x - x_0) ** 2
frac = - xdiff / (4 * (sigma ** 2))
psi_init = np.exp(frac + 1.0j * kappa * x)

# Compute psi0
psi_squared = psi_init * np.conj(psi_init)
psi0_squared = 1 / (np.sum(psi_squared) * a)
psi0 = np.sqrt(psi0_squared)

# Compute psi initially
psi = psi_init * psi0


# Define square well potential
def V_sqaure_well(*args):
    return 0


# Define Quantum Harmonic Oscillator potential
def V_QHO(x, m=m, w=3e15):
    return 0.5 * m * (w ** 2) * (x ** 2)


# Define double well potential
def V_double_well(x, x_1=L / 3, V_0=6e-17):
    frac = (x ** 2) / (x_1 ** 2)
    return V_0 * (frac - 1) ** 2


# Define H_D() Matrix
def H_D(p=p, P=P, L=L):
    A = -(hbar ** 2) / (2 * m * (a ** 2))
    B_p = V_sqaure_well(p * a - L / 2) - 2 * A
    diagonal = np.eye(P - 1, k=0) * B_p
    sup_diag = np.eye(P - 1, k=1) * A
    sub_diag = np.eye(P - 1, k=-1) * A
    return diagonal + sub_diag + sup_diag


# Define energy expectation
def energyExpectation(psi):
    psi_conj = np.conjugate(psi)
    H_times_ket = np.dot(H_D(), psi)
    bra = psi_conj
    return np.sum(np.dot(bra, H_times_ket)) * a


# Set initial time step
tstep = 0

# Set initial value of psi
psivec = np.array(psi[:P - 1], complex)

# Set expectation of position data
xdata = []

# Set normalization data
normdata = []

# Set energy expectation data
energydata = []

# Set total time to iterate
T = N * tau

# Define Identity Matrix
I = np.identity(P - 1, complex)

# Compute Lmatrix and Rmatrix
H = H_D()
Lmatrix = I + 1.0j * (tau / (2 * hbar)) * H
Rmatrix = I - 1.0j * (tau / (2 * hbar)) * H

# Set current psi to initial psi
psicurr = psivec

# Crank-Nicolson loop
while tstep < T:
    print(f"[STATUS] Solving PDE... {round(100 * tstep / T, 1)}%")

    # compute v and solve for new psi
    vvec = np.matmul(Rmatrix, psicurr)
    psinew = np.linalg.solve(Lmatrix, vvec)

    # Append values of expectations as required.
    xdata.append(abs(positionExpectation(psinew, x[:P - 1])))
    energydata.append(abs(energyExpectation(psinew)))
    normdata.append(np.sum(abs(psi) ** 2) * a)

    # Save psi's for T, T/2, T/4
    if np.isclose(tstep, T / 4, atol=tau):
        psi_quater = psinew

    if np.isclose(tstep, T / 2, atol=tau):
        psi_half = psinew

    if np.isclose(tstep, T, atol=tau):
        psi_T = psinew

    # Update time step
    tstep += tau

    # Update psi for next loop
    psicurr = psinew

print("\n[STATUS] Solved!")
print("[STATUS] Plotting figures now...")

# Set time array for plotting purposes
t = np.arange(0, T, tau)

part = "a"

# Plot real value of psi
plt.figure(figsize=(8, 6))
plt.plot(x[:P - 1], np.real(psivec))
plt.title("Normalized wavefunction with Respect to Position at $t=0$")
plt.xlabel("$x(p) = a * p$")
plt.ylabel("Intensity")
plt.savefig(f"Q1{part}-1.pdf")

# Plot all psi's evolved in time.
plt.figure(figsize=(8, 6))
plt.plot(x[:P - 1], abs(psivec * np.conj(psivec)), label='t = 0')
plt.plot(x[:P - 1], abs(psi_half) ** 2, label="t = T/2")
plt.plot(x[:P - 1], abs(psi_quater) ** 2, label="t = T/4")
plt.plot(x[:P - 1], abs(psi_T) ** 2, label="t = T")
plt.title("Normalized wavefunction with Respect to Position")
plt.xlabel("$x(p) = a * p$")
plt.ylabel("Intensity")
plt.legend()
plt.savefig(f"Q1{part}-2.pdf")
plt.close()

# Plot values of energy expectation
energydata.pop()
plt.figure()
plt.plot(t, energydata)
plt.title("Energy with respect to time")
plt.xlabel("Time (t)")
plt.ylabel("Energy")
plt.ylim(-1e-16, 1e-16)
plt.savefig(f"Q1{part}-3.pdf")

# Plot values of position expectation
xdata.pop()
plt.figure()
plt.plot(t, xdata)
plt.title("Position Expectation with respect to time")
plt.xlabel("Time (t)")
plt.ylabel("Position")
plt.savefig(f"Q1{part}-4.pdf")
plt.ylim(min(xdata), max(xdata))

# Plot values of normed data
normdata.pop()
plt.figure()
plt.plot(t, normdata)
plt.title("Normalization of Psi with respect to time")
plt.xlabel("Time (t)")
plt.ylabel("Normalization")
plt.savefig(f"Q1{part}-5.pdf")

# Be happy
print("[STATUS] All Done!")
