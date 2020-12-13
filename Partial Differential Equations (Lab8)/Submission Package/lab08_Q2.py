# Author: Utkarsh MAli
import numpy as np
import matplotlib.pyplot as plt

# Question 2b
# Set initial conditions
g = 9.81
J = 51
dx = 0.02
L = 1.
x = np.arange(0, L + dx, dx)
n_b = 0.
H = 0.01
dt = 0.01
A = 0.002
mu = 0.5
sigma = 0.05
tf = 4.5
t = 0.

# Initialize empty u, u+1 and u-1 arrays
u_new = np.zeros(x.shape)
u = np.zeros(x.shape)

# Same array initialization for n
num = -(x - mu) ** 2
dem = sigma ** 2
factor = num / dem
exp_power = A * np.exp(factor)
nstart = H + exp_power - np.mean(exp_power)
eta_new = np.zeros(x.shape)
eta = nstart


# Define F equation (7), assume n_b = 0 be default.
def F(u, n, g=9.81, n_b=0.):
    return 0.5 * (u ** 2) + g * n, (n - n_b) * u


def boundary_condition(u, condition=0):
    if np.isclose(u[0], condition):
        u[0] = condition
    if np.isclose(u[-1], condition):
        u[-1] = condition
    return u


# Initialize figure
plt.figure(figsize=(8, 6))
while t < tf:
    newF = F(u, eta)

    # Handle edge cases first with forward/backward difference.
    eta_new[0] = eta[0] - (dt / dx) * (newF[1][1] - newF[1][0])
    eta_new[J - 1] = eta[J - 1] - (dt / dx) * (newF[1][J - 1] - newF[1][J - 2])

    # Cover the rest in a for loop
    for j in range(1, J - 1):
        eta_new[j] = eta[j] - (dt / (2 * dx)) * (newF[1][j + 1] - newF[1][j - 1])
        u_new[j] = u[j] - (dt / (2 * dx)) * (newF[0][j + 1] - newF[0][j - 1])

    u_new = boundary_condition(u_new)
    u = np.copy(u_new)
    eta = np.copy(eta_new)

    for t_wanted in [0, 1, 4]:
        if np.isclose(t, t_wanted):
            # if abs(t-t_wanted) < 1e-10:
            plt.plot(x, eta, label=f"$t =${round(t, 2)}s")

    t += dt

plt.xlabel("Position $x$")
plt.xlim(0, 1)
plt.title("Implementing the 1D shallow water system with the FTCS scheme")
plt.ylabel("Altitiude of free surface $\eta(x,t)$")
plt.legend()
plt.savefig("Q2a.pdf")
plt.show()
