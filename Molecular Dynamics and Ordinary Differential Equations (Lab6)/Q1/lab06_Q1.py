# Author: Utkarsh Mali
import numpy as np
import matplotlib.pyplot as plt

# Question 1b
""" Using RK4 to Solve
"""

# Define units
G = 1.
M = 10.
L = 2.


# Define functions, after reducing them to first order
def f(r, t):
    """
    This function computes 2 second order linear differential equations by collapsing them into 4 simpler first order
    simple differential equations. We use the terms xdot and ydot to represent 2 of our collapsed simple differential
    equations.
    :param r: Distance and Speed vector
    :param t: Time
    :return: Array of r vector
    """
    # We first set our r vector equal to our position and velocities as required
    x, xdot, y, ydot = r[0], r[1], r[2], r[3]

    # We then calculate our distance r as defined in the question
    d = np.sqrt(x ** 2 + y ** 2)

    # Then set parameters of constants to be multiplied to the numerator
    c = - G * M

    # Depending on which equation is being solved, set the value of the numerator to xc and yc.
    num1 = x * c
    num2 = y * c

    # Set the appropriate value of denominator as defined in the question
    dem = d ** 2 * np.sqrt(d ** 2 + (L ** 2) / 4)

    # compute fraction for x ODE
    x_eq = num1 / dem

    # compute fraction for y ODE
    y_eq = num2 / dem

    # return array of updated positions and velocities.
    return np.array([xdot, x_eq, ydot, y_eq])


# Set initial conditions of time
t0 = 0.0
tend = 10.0

# Set initial conditions of position and velocity
x0, y0 = 1.0, 0.0
xdot, ydot = 0.0, 1.0

# Set N and h
N = 10000
h = (tend - t0) / N

# Setting arrays for iteration
tpoints = np.arange(t0, tend, h)
xpoints = []
xdotpoints = []
ypoints = []
ydotpoints = []

# Set r vector
r = np.array([x0, xdot, y0, ydot])

# Apply RK4 method
for t in tpoints:
    # append points first
    xpoints.append(r[0])
    xdotpoints.append(r[1])
    ypoints.append(r[2])
    ydotpoints.append(r[3])

    # Compute RK4 method loop
    # Notice all the t's being used are dummy variables which do not affect the computation
    k1 = h * f(r, t)
    k2 = h * f(r + 0.5 * k1, t + 0.5 * h)
    k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
    k4 = h * f(r + k3, t + h)

    # update r for next iteration
    r += (k1 + 2 * k2 + 2 * k3 + k4) / 6

# Plotting the figure requested in the question
plt.figure()
plt.plot(xpoints, ypoints, color='firebrick')
plt.xlabel("Position $x$")
plt.ylabel("Position $y$")
plt.title("Position of Orbiting Space Garbage")
plt.savefig("fig1.pdf")

# Plot an additional speed figure for aesthetic
plt.figure()
plt.plot(tpoints, xdotpoints, label="Velocity $x$", linestyle = '--')
plt.plot(tpoints, ydotpoints, label="Velocity $y$", linestyle = ':', color = 'green')
plt.plot(tpoints, np.sqrt(np.array(xdotpoints) ** 2 + np.array(ydotpoints) ** 2), label="Velocity $r$", color='salmon')
plt.legend()
plt.title("Speed as a function of time of Orbiting Space Garbage")
plt.xlabel("Time $s$")
plt.ylabel("Speed $m/s$")
plt.savefig('fig2.pdf')
