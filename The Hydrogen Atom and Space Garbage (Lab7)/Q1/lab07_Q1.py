# Author: Utkarsh Mali

import numpy as np
import time
import matplotlib.pyplot as plt

# Define units
G = 1.
M = 10.
L = 2.
print("[STATUS] Initializing...")


# Copying defined function taken from Lab 6, Question 1
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


# Question 1a
# Time the run with uniform time intervals
uniform_t0 = time.time()
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
plt.plot(xpoints, ypoints, color='firebrick', label='Uniform Interval')

# Calculate final time
uniform_tf = time.time() - uniform_t0

# Print time taken
print("[STATUS] Running...")
print("Uniform Interval Runtime:", round(uniform_tf, 4), "sec")

# Question 1b
# Time the varying time step method.
varying_t0 = time.time()

# Set initial parameters, using same ones as before.
delta = 10e-6
r = np.array([x0, xdot, y0, ydot])

# Set N and h
N = 10000
h = (tend - t0) / N

# Setting arrays for iteration
tpoints = []
xpoints = []
xdotpoints = []
ypoints = []
ydotpoints = []


# Define rho function
def rho(h, delta, x1, x2, y1, y2):
    # Calculate epsilon
    epsilon_x = (1 / 30) * (x1 - x2)
    epsilon_y = (1 / 30) * (y1 - y2)

    # Multiply by delta
    num = h * delta
    temp1 = epsilon_x ** 2
    temp2 = epsilon_y ** 2

    # Calculate Norm
    dem = np.sqrt(temp1 + temp2)
    return num / dem


# Set initial time
t = t0

# Apply RK4 method
while t < tend:
    # First do 2 steps of the solution, each of size h
    # Doing first step
    k1 = h * f(r, t)
    k2 = h * f(r + 0.5 * k1, t + 0.5 * h)
    k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
    k4 = h * f(r + k3, t + h)
    rtemp0 = r + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    # Set new t and complete first step again
    tnew = t + h
    k1 = h * f(rtemp0, tnew)
    k2 = h * f(rtemp0 + 0.5 * k1, tnew + 0.5 * h)
    k3 = h * f(rtemp0 + 0.5 * k2, tnew + 0.5 * h)
    k4 = h * f(rtemp0 + k3, tnew + h)
    rtemp1 = rtemp0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    # Runge-Kutta step but with 2 * h
    hnew = 2 * h
    k1 = hnew * f(r, t)
    k2 = hnew * f(r + 0.5 * k1, t + 0.5 * hnew)
    k3 = hnew * f(r + 0.5 * k2, t + 0.5 * h)
    k4 = hnew * f(r + k3, t + hnew)
    rtemp2 = r + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    # Compute rho
    rho_ = rho(h, delta, rtemp1[0], rtemp2[0], rtemp1[2], rtemp2[2])

    # Compute new h
    h *= rho_ ** (1 / 4)

    # Limiting constant of h as defined in Newman
    h = min(h, 0.015)

    # Add values for large rho, repeat calculation for smaller rho
    if rho_ > 1.:

        # Append time to time list
        tpoints.append(t)

        # Change time according to computed varying h
        t += 2 * h

        # Set r
        r = rtemp1

        # Update values of r
        xpoints.append(r[0])
        xdotpoints.append(r[1])
        ypoints.append(r[2])
        ydotpoints.append(r[3])

# Print time taken
varying_tf = time.time() - varying_t0
print("Varying Interval Runtime:", round(varying_tf, 4), "sec")

# Plotting the figure requested in the question
plt.plot(xpoints, ypoints, 'k.', label='Varying Interval')
plt.xlabel("Position $x$")
plt.ylabel("Position $y$")
plt.title("Position of Orbiting Space Garbage")
plt.legend()
plt.savefig('fig1.pdf')

# Question 1c
plt.figure()
dtpoints = np.array(tpoints[1:]) - np.array(tpoints[:-1])
plt.plot(tpoints[:-1], dtpoints)  # drop the last point in tpoints
plt.xlabel("Time ($s$)")
plt.ylabel("Time Step ($\Delta t$)")
plt.title("Size of time steps as a function of time")
plt.savefig("fig2.pdf")

# Plot an additional speed figure for aesthetic
plt.figure()
plt.plot(tpoints, xdotpoints, label="Velocity $x$", linestyle='--')
plt.plot(tpoints, ydotpoints, label="Velocity $y$", linestyle=':', color='green')
plt.plot(tpoints, np.sqrt(np.array(xdotpoints) ** 2 + np.array(ydotpoints) ** 2), label="Velocity $r$", color='salmon')
plt.legend()
plt.title("Speed as a function of time of Orbiting Space Garbage")
plt.xlabel("Time $s$")
plt.ylabel("Speed $m/s$")
plt.savefig('fig3.pdf')

print("\n \n[STATUS] All done!")
