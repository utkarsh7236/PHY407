# Author: Utkarsh Mali

import numpy as np
import matplotlib.pyplot as plt
import random as rn

# Ex 10.3: Brownian Motion
print('[STATUS] Initializing...')

# Define Constants
L = 101
steps = 5000

# Set i, j
i = np.arange(0, L, 1)
j = np.arange(0, L, 1)

# Define starting position
x0, y0 = 50, 50

# Define empty arrays for x and y.
xdata = []
ydata = []


# Use nextmove defenition from Prof Nicolas
def nextmove(x, y):
    """ randomly choose a direction
    0 = up, 1 = down, 2 = left, 3 = right"""
    direction = rn.randrange(0, 4)

    if direction == 0:  # move up
        y += 1
    elif direction == 1:  # move down
        y -= 1
    elif direction == 2:  # move right
        x += 1
    elif direction == 3:  # move left
        x -= 1
    else:
        print("[ERROR] Direction isn't 0-3")

    return x, y


# Define boundary checker
def check_boundary(x, y, L=L):
    retx, rety = x, y
    if x > (L - 1):
        retx = L - 1
    elif x < 0:
        retx = 0
    if y > (L - 1):
        rety = L - 1
    elif y < 0:
        rety = 0
    return retx, rety

# Set current x and y in loop to initial.
curr_x = x0
curr_y = y0

# Run simulation through for loop
print("[STATUS] Running Monte-Carlo Simulation...")
for index in range(steps):
    # set temp x and y with nextmove()
    temp_x, temp_y = nextmove(curr_x, curr_y)

    # set new x and y after checking boundary.
    new_x, new_y = check_boundary(temp_x, temp_y)

    # Update xdata and ydata
    xdata.append(curr_x)
    ydata.append(curr_y)

    # update current x and y in loop
    curr_x, curr_y = new_x, new_y

# Plot figure, set to square and limit axis
plt.figure(figsize=(8, 8))
plt.title("Complete Random Walk of a Particle (5000 steps)")
plt.plot(xdata, ydata)
plt.xlim(0, max(i))
plt.ylim(0, max(j))
plt.xlabel("Position $x$")
plt.ylabel("Position $y$")
plt.savefig("q1a.pdf")

print("[STATUS] All done!")
