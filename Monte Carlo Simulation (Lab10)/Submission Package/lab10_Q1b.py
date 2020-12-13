# Author: Utkarsh Mali
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
import random as rn

# Ex 10.13: Diffusion-limited aggregation
print('[STATUS] Initializing...')

# Define Constants
L = 101
steps = 5000

# Set i, j
i = np.arange(0, L, 1)
j = np.arange(0, L, 1)

# Define starting position
x0, y0 = (L - 1) / 2, (L - 1) / 2


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
def atBoundary(x, y, L=L):
    if x >= L:
        return True
    elif x <= 0:
        return True
    if y >= L:
        return True
    elif y <= 0:
        return True
    return False


# Define anchor checker
def atAnchor(x, y, anchorz):
    for anchor in anchorz:
        if x == anchor[0] and y == anchor[1]:
            return True
        else:
            pass
    return False


# Run simulation through for loop
print("[STATUS] Running Monte-Carlo Simulation...")

# Set number of particles
numParticles = 100

# Set list for anchored points
anchors = []

# For loop over all the particles
for k in range(numParticles + 1):
    # Print status update every few loops.
    if k % 20 == 0:
        print(f"[STATUS] Monte-Carlo Progress: {100 * k / numParticles}%")

    # Empty xdata and ydata for each particle
    xdata = []
    ydata = []

    # Set each particle to the center at the start
    curr_x = x0
    curr_y = y0

    # Run Monte-Carlo Simulation for single particle
    # Loop over the number of steps we want each particle to take
    for index in range(steps):

        # Find new values of x and y.
        new_x, new_y = nextmove(curr_x, curr_y)

        # Append to xdata and ydata
        xdata.append(curr_x)
        ydata.append(curr_y)

        # If item is at anchor then we can stop looping
        if atAnchor(new_x, new_y, anchors):
            # Add location of particce to all anchors
            anchors.append([curr_x, curr_y])
            break

        # Update current x and y and continue simualtion
        curr_x, curr_y = new_x, new_y

        # If item is at boundary, we can stop looping
        if atBoundary(new_x, new_y):
            # Add location of particle to all anchors
            anchors.append([new_x, new_y])
            break

print("\n[STATUS] Simulated!")
print("[STATUS] Plotting Figures...")
print("[STATUS} No figures to plot! Cleaning up...")
print("[STATUS] All done!")
