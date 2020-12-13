# Header
# Author Utkarsh Mali

import numpy as np
import matplotlib.pyplot as plt

print("[STATUS] COMPILING")

# Pseudo Code Q1b #
""" Write a “pseudo-code” for a program that integrates (using Euler-Cromer) your 
equations to calculate the position and velocity of the planet as a function of 
time under the Newtonian gravity force. The output of your code should include 
graphs of the components of velocity as a function of time and a plot of the 
orbit (x vs. y) in space.
"""
# gravitational constant (in terms of solar mass and au)
# initial x position
# initial y position
# initial x velocity
# initial y velocity
# time step
# empty array for x position
# empty array for y position
# empty array for time after time step
# empty array for x velocity
# empty array for y velocity
# empty array for angular momentum
# add initial x position to x array
# add initial y position to y array
# add initial time 0 to time array
# add initial x velocity to vx array
# add initial y velocity to vy array
# while loop calculating time stepped x and y positions and velocities (loop keeps running while time <= 1 year)
#       calculate and update angular momentum
#       update vx using equation 6a
#       add updated vx to vx array
#       use updated vx to update x (Euler-Cromer method)
#       update vy using equation 6b
#       add updated vy to vy array
#       use updated vy to update y (Euler-Cromer method)
#       update x,y arrays with new values
#       update time with time step
#       add updated time to time array
# new figure
# plot appended x and y array in 2d space
# label axis
# new figure
# plot appended vx and vx as a function of time
# label axis

# Code to Question 1c #
"""Use a time step ∆t = 0.0001 yr and integrate for 1 year. Check if angular 
momentum is conserved from the beginning to the end of the orbit. You should 
see an elliptical orbit.
"""
print("[STATUS] EXECUTING 1C")
# gravitational constant (in terms of solar mass and au)
G = 39.5

# initial x position
x = 0.47

# initial y position
y = 0.0

# initial x velocity
vx = 0.0

# initial y velocity
vy = 8.17

# time step
dt = 0.0001

# empty array for x position, add initial x position to x array
x_pos = np.array([x])

# empty array for y position, add initial y position to y array
y_pos = np.array([y])

# empty array for time after time step, add initial time 0 to time array
t = 0
time = np.array([t])

# empty array for x velocity, add initial x velocity to vx array
vx_arr = np.array([vx])

# empty array for y velocity, add initial y velocity to vy array
vy_arr = np.array([vy])

# empty arr for angular momentum
merc_mass = 1.651 * 10**(-7)
L_arr = np.array([merc_mass * np.sqrt(vx**2 + vy ** 2) * np.sqrt(x ** 2 + y ** 2)])

# loop update x and y positions and velocities (while time <= 1 year)
while t <= 1:

    L = merc_mass * np.sqrt(vx**2 + vy ** 2) * np.sqrt(x ** 2 + y ** 2)
    L_arr = np.append(L_arr, L)

    # update vx using equation 6a
    vx = vx - (G * 1 * x * dt) / (np.sqrt(x ** 2 + y ** 2)) ** 3

    # add updated vx to vx array
    vx_arr = np.append(vx_arr, vx)

    # update vy using equation 6b
    vy = vy - (G * 1 * y * dt) / (np.sqrt(x ** 2 + y ** 2)) ** 3

    # add updated vy to vy array
    vy_arr = np.append(vy_arr, vy)

    # use updated vy to update y (Euler-Cromer method)
    x = x + vx * dt

    # use updated vx to update x (Euler-Cromer method)
    y = y + vy * dt

    # update x,y arrays with new values
    x_pos = np.append(x_pos, x)
    y_pos = np.append(y_pos, y)

    # update time with time step
    t = t + dt

    # add updated time to time array
    time = np.append(time, t)

# print(vx_arr)
# print(vy_arr)
# print(time)
# print(x_pos)
# print(y_pos)

# new figure
plt.figure()

# plot appended x and y array in 2d space
plt.plot(x_pos, y_pos)

# label axis
scale = 1.1
xmin, xmax = -max(x_pos)*scale, max(x_pos)*scale
ymin, ymax = -max(y_pos)*scale, max(y_pos)*scale
plt.title("Mercury orbit for 1 year plotted in space")
plt.xlabel("Position x (AU)")
plt.ylabel("Position y (AU)")
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.grid()
plt.savefig("fig1.pdf")

# new figure
plt.figure()

# plot appended vx and vx as a function of time
plt.plot(time, vx_arr, label="Velocity x")
plt.plot(time, vy_arr, label="Velocity y")

# label axis
plt.legend()
plt.title("Bi-directional velocity as a function of time for 1 year")
plt.xlabel("Time (yr)")
plt.ylabel("Magnitude of Velocity (AU/yr)")
plt.grid()
plt.savefig("fig2.pdf")

print(f"Angular momentum was always between {round(max(L_arr),10)} and {round(min(L_arr), 10)}"
      f" hence not fully conserved. \nWe do not expect the angular momentum in"
      f" eliptical orbits to"
      f" always be conserved since \nthere is a sin factor between the angle to"
      f" the center of the orbit and the actual \nperpendicular vector to v.")

# Question 1d #
"""Now alter the gravitational force in your code to the general relativity form given in
eqn. (7). The actual value of α for Mercury is given in the physics background, but
it will be too small for our computational framework here. Instead, try α = 0.01 AU2
which will exaggerate the effect. Demonstrate Mercury’s orbital precession by plotting
several orbits in the x, y plane that show the perihelion (furthest point) of the orbit
moves around in time.
"""
print("[STATUS] EXECUTING 1D")
# gravitational constant (in terms of solar mass and au)
G = 39.5

# initial x position
x = 0.47

# initial y position
y = 0.0

# initial x velocity
vx = 0.0

# initial y velocity
vy = 8.17

# time step
dt = 0.0001

# alpha
a = 0.01

# empty array for x position, add initial x position to x array
x_pos = np.array([x])

# empty array for y position, add initial y position to y array
y_pos = np.array([y])

# empty array for time after time step, add initial time 0 to time array
t = 0
time = np.array([t])

# empty array for x velocity, add initial x velocity to vx array
vx_arr = np.array([vx])

# empty array for y velocity, add initial y velocity to vy array
vy_arr = np.array([vy])

# empty arr for angular momentum
merc_mass = 1.651 * 10**(-7)
L_arr = np.array([merc_mass * np.sqrt(vx**2 + vy ** 2) * np.sqrt(x ** 2 + y ** 2)])

# loop update x and y positions and velocities (while time <= 1 year)
while t <= 2:

    L = merc_mass * np.sqrt(vx**2 + vy ** 2) * np.sqrt(x ** 2 + y ** 2)
    L_arr = np.append(L_arr, L)

    # update vx using equation 6a
    r = np.sqrt(x ** 2 + y ** 2)
    vx = vx - ((G * 1 * x * dt) / (r ** 3)) * (1 + a/(r**2))

    # add updated vx to vx array
    vx_arr = np.append(vx_arr, vx)

    # update vy using equation 6b
    vy = vy - ((G * 1 * y * dt) / (r ** 3)) * (1 + a/(r**2))

    # add updated vy to vy array
    vy_arr = np.append(vy_arr, vy)

    # use updated vy to update y (Euler-Cromer method)
    x = x + vx * dt

    # use updated vx to update x (Euler-Cromer method)
    y = y + vy * dt

    # update x,y arrays with new values
    x_pos = np.append(x_pos, x)
    y_pos = np.append(y_pos, y)

    # update time with time step
    t = t + dt

    # add updated time to time array
    time = np.append(time, t)

# new figure
plt.figure()

# plot appended x and y array in 2d space
plt.plot(x_pos, y_pos)

# label axis
scale = 1.1
xmin, xmax = -max(x_pos)*scale, max(x_pos)*scale
ymin, ymax = -max(y_pos)*scale, max(y_pos)*scale
plt.title("GR Mercury orbit for 2 years plotted in space")
plt.xlabel("Position x (AU)")
plt.ylabel("Position y (AU)")
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.grid()
plt.savefig("fig3.pdf")

# new figure
plt.figure()

# plot appended vx and vx as a function of time
plt.plot(time, vx_arr, label="Velocity x")
plt.plot(time, vy_arr, label="Velocity y")

# label axis
plt.legend()
plt.title("GR Bi-directional velocity as a function of time for 1 year")
plt.xlabel("Time (yr)")
plt.ylabel("Magnitude of Velocity (AU/yr)")
plt.grid()
plt.savefig("fig4.pdf")

print("[STATUS] FINISHED")
