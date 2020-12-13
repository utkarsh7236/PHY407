# Author: Utkarsh Mali
import matplotlib.pyplot as plt
import numpy as np

# QUESTION 2A
# Setting constants, some are taken from Newman812a-sol.py provided by Prof Grisouard
# [kg] Mass of the Sun
M = 1.9891e30

# Constant to convert to year time
sec_to_year_convert = 24 * 365.25 * 60 * 60

# Gravitational Constant
G = 6.6738e-11 * sec_to_year_convert ** 2

# Divide by 53 since we are working in years
H = 1 / 53

# Iniital x
x0 = 0

# Initial y
y0 = 1.4710e11

# Initial xdot
xdot0 = 3.0287e4 * sec_to_year_convert

# Initial ydot
ydot0 = 0

# Accuracy given in the question
delta = 1e4


# Define Function
def f(r):
    x, xdot = r[0], r[1]
    y, ydot = r[2], r[3]

    def xpos(x, y):
        return -G * M * x / (np.sqrt(x ** 2 + y ** 2) ** 3)

    def ypos(x, y):
        return -G * M * y / (np.sqrt(x ** 2 + y ** 2) ** 3)

    return np.array([xdot, xpos(x, y), ydot, ypos(x, y)])


# Set empty arrays
xpoints, ypoints = [], []
xdotpoints, ydotpoints = [], []

# Set initializing vector
r = np.array([x0, xdot0, y0, ydot0])

# Do the "big steps" of size H, looping for 1 year
# Loop over using Newman's Bulirsch-Stoer Method
for t in range(54):

    # Update initial points into final arrays
    xpoints.append(r[0])
    xdotpoints.append(r[1])
    ypoints.append(r[2])
    ydotpoints.append(r[3])

    # Do one modified midpoint step to get things started
    n = 1
    r1 = r + 0.5 * H * f(r)
    r2 = r + H * f(r1)

    # The array R1 stores the first row of the
    # extrapolation table, which contains only the single
    # modified midpoint estimate of the solution at the
    # end of the interval
    # Set the empty dimension to 4
    R1 = np.empty([1, 4], float)
    R1[0] = 0.5 * (r1 + r2 + 0.5 * H * f(r2))

    # Now increase n until the required accuracy is reached
    error = 2 * H * delta
    while error > H * delta:

        n += 1
        h = H / n

        # Modified midpoint method
        r1 = r + 0.5 * h * f(r)
        r2 = r + h * f(r1)
        for i in range(n - 1):
            r1 += h * f(r2)
            r2 += h * f(r1)

        # Calculate extrapolation estimates.  Arrays R1 and R2
        # hold the two most recent lines of the table
        R2 = R1

        # Set the empty dimension to 4
        R1 = np.empty([n, 4], float)
        R1[0] = 0.5 * (r1 + r2 + 0.5 * h * f(r2))
        for m in range(1, n):
            epsilon = (R1[m - 1] - R2[m - 1]) / ((n / (n - 1)) ** (2 * m) - 1)
            R1[m] = R1[m - 1] + epsilon

        # Take the error of both x and y values.
        error = abs(np.sqrt(epsilon[0] ** 2 + epsilon[2] ** 2))

    # Set r equal to the most accurate estimate we have,
    # before moving on to the next big step
    r = R1[n - 1]

# Plot the results
plt.plot(xpoints, ypoints, color='darkorange')
v = 1.6e11
plt.xlim(-v, v)
plt.ylim(-v, v)
plt.xlabel("Position $x$")
plt.ylabel("Position $y$")
plt.title("Earths Orbit using Bulirsch-Stoer")
plt.savefig("fig2a.pdf")

# QUESTION 2B
# Setting constants, like previously. Change values of H, and initial setup
H = 0.5
x0 = 4.4368e12
y0 = 0
xdot0 = 0
ydot0 = 6.1218e3 * sec_to_year_convert

# Set empty arrays for Pluto's Orbit
xpoints, ypoints = [], []
xdotpoints, ydotpoints = [], []

# Initialize vector
r = np.array([x0, xdot0, y0, ydot0])

# Loop over using Newman's Bulirsch-Stoer Method
for t in range(500):

    # Set initial values of x, y, xdot, ydot
    xpoints.append(r[0])
    xdotpoints.append(r[1])
    ypoints.append(r[2])
    ydotpoints.append(r[3])

    # Do one modified midpoint step to get things started
    n = 1
    r1 = r + 0.5 * H * f(r)
    r2 = r + H * f(r1)

    # The array R1 stores the first row of the
    # extrapolation table, which contains only the single
    # modified midpoint estimate of the solution at the
    # end of the interval
    R1 = np.empty([1, 4], float)
    R1[0] = 0.5 * (r1 + r2 + 0.5 * H * f(r2))

    # Now increase n until the required accuracy is reached
    error = 2 * H * delta
    while error > H * delta:

        n += 1
        h = H / n

        # Modified midpoint method
        r1 = r + 0.5 * h * f(r)
        r2 = r + h * f(r1)
        for i in range(n - 1):
            r1 += h * f(r2)
            r2 += h * f(r1)

        # Calculate extrapolation estimates.  Arrays R1 and R2
        # hold the two most recent lines of the table
        R2 = R1
        R1 = np.empty([n, 4], float)
        R1[0] = 0.5 * (r1 + r2 + 0.5 * h * f(r2))
        for m in range(1, n):
            epsilon = (R1[m - 1] - R2[m - 1]) / ((n / (n - 1)) ** (2 * m) - 1)
            R1[m] = R1[m - 1] + epsilon
        error = abs(np.sqrt(epsilon[0] ** 2 + epsilon[2] ** 2))

    # Set r equal to the most accurate estimate we have,
    # before moving on to the next big step
    r = R1[n - 1]

# Plot the results
plt.figure()
plt.plot(xpoints, ypoints, color='firebrick')
plt.title("Pluto's Orbit using Bulirsch-Stoer")
plt.xlabel("Position $x$")
plt.ylabel("Position $y$")
v = 8e12
plt.xlim(-v, v)
plt.ylim(-v, v)
plt.savefig("fig2b.pdf")
