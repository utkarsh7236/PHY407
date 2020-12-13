import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# for formatting the graph
sns.set()

# set constants
sigma = 1
epsilon = 1
mass = 1.0 # mass of the molecule


# x and y components of accelerations found in part (a)
def acceleration_x(r1x, r1y, r2x, r2y):
    """
    Return aceleration array with components in x direction and y direction.
    r1, r2: lists with two elements i.e. coordinates (x,y)
    """
    # the distance between the two particles
    r = np.sqrt((r2y - r1y)**2 + (r2x - r1x)**2)
    # use F = V(r) = ma where m = 1 to find a. 
    # Multiply by x2-x1 / r for x component
    acceleration_x = 24 * (r2x - r1x) * (-2 * r**(-13) + r**(-7))
    return acceleration_x

def acceleration_y(r1x, r1y, r2x, r2y):
    """
    Return aceleration array with components in x direction and y direction.
    r1, r2: lists with two elements i.e. coordinates (x,y)
    """
    # the distance between the two particles
    r = np.sqrt((r2y - r1y)**2 + (r2x - r1x)**2)
    # use F = V(r) = ma where m = 1 to find a. 
    # Multiply by y2-y1 / r for y component
    acceleration_y = 24 * (r2y - r1y) * (-2 * r**(-13) + r**(-7))
    return acceleration_y


def run_Verlet(initial_r1, initial_r2):
    """
    r1, r2: lists with two elements i.e. coordinates (x,y)
    """
    
    # set given values
    h = 0.01 # dt
    num_steps = 100
    
    # initialize arrays to hold results
    r1x = np.zeros(num_steps)
    r1y = np.zeros(num_steps)
    r2x = np.zeros(num_steps)
    r2y = np.zeros(num_steps)
    v1x = np.zeros(num_steps)
    v1y = np.zeros(num_steps)
    v2x = np.zeros(num_steps)
    v2y = np.zeros(num_steps)
    
    r1x[0] = initial_r1[0]
    r1y[0] = initial_r1[1]
    r2x[0] = initial_r2[0]
    r2y[0] = initial_r2[1]
    
    # for the first step only. Represent v(t + h/2) in eq7 by v_temp
    v1x_temp = v1x[0] + 1/2 * h * acceleration_x(initial_r1[0], initial_r1[1],
                  initial_r2[0], initial_r2[1])
    v1y_temp = v1y[0] + 1/2 * h * acceleration_y(initial_r1[0], initial_r1[1],
                  initial_r2[0], initial_r2[1])
    v2x_temp = v2x[0] + 1/2 * h * -acceleration_x(initial_r1[0], initial_r1[1],
                  initial_r2[0], initial_r2[1])
    v2y_temp = v2y[0] + 1/2 * h * -acceleration_y(initial_r1[0], initial_r1[1],
                  initial_r2[0], initial_r2[1])
    
    # apply eqn 8 - 11 repeatedly
    for i in range(1, num_steps):
        # apply equation 8
        r1x[i] = r1x[i-1] + h * v1x_temp
        r1y[i] = r1y[i-1] + h * v1y_temp
        r2x[i] = r2x[i-1] + h * v2x_temp
        r2y[i] = r2y[i-1] + h * v2y_temp
        
        # apply equation 9. K is a list [kx, ky]
        k1x = h * acceleration_x(r1x[i], r1y[i], r2x[i], r2y[i]) 
        k1y = h * acceleration_y(r1x[i], r1y[i], r2x[i], r2y[i])
        k2x = h * -acceleration_x(r1x[i], r1y[i], r2x[i], r2y[i])
        k2y = h * -acceleration_y(r1x[i], r1y[i], r2x[i], r2y[i])
        
        # apply equation 10
        v1x[i] = v1x_temp + 1/2 * k1x
        v1y[i] = v1y_temp + 1/2 * k1y
        v2x[i] = v2x_temp + 1/2 * k2x
        v2y[i] = v2y_temp + 1/2 * k2y
           
        # apply equation 11
        v1x_temp = v1x_temp + k1x
        v1y_temp = v1y_temp + k1y
        v2x_temp = v2x_temp + k2x
        v2y_temp = v2y_temp + k2y
        
    return r1x, r1y, r2x, r2y, v1x, v1y, v2x, v2y
    

# First initial condition q (i)
r1_i, r2_i = [4, 4], [5.2, 4]
r1x_i, r1y_i, r2x_i, r2y_i, v1x_i, v1y_i, v2x_i, v2y_i = run_Verlet(r1_i, r2_i) 
plt.figure(figsize=(8,6))
plt.plot(r1x_i, r1y_i, '.', label="particle 1")
plt.plot(r2x_i, r2y_i, '.', label="particle 2" )
plt.legend(fontsize=14)
plt.xlabel("x(t)", fontsize=14)
plt.ylabel("y(t)", fontsize=14)
plt.title("Trajectories of particles for the initial condition (i)", fontsize=16)
plt.savefig("2a.pdf")


# Second initial condition q (ii)
r1_ii, r2_ii = [4.5, 4], [5.2, 4]
r1x_ii, r1y_ii, r2x_ii, r2y_ii, v1x_ii, v1y_ii, v2x_ii, v2y_ii = run_Verlet(r1_ii, r2_ii) 
plt.figure(figsize=(8,6))
plt.plot(r1x_ii, r1y_ii, '.', label="particle 1")
plt.plot(r2x_ii, r2y_ii, '.', label="particle 2" )
plt.legend(fontsize=14)
plt.xlabel("x(t)", fontsize=14)
plt.ylabel("y(t)", fontsize=14)
plt.title("Trajectories of particles for the initial condition (ii)", fontsize=16)
plt.savefig("2b.pdf")


# Third initial condition q (iii)
r1_iii, r2_iii = [2, 3], [3.5, 4.4]
r1x_iii, r1y_iii, r2x_iii, r2y_iii, v1x_iii, v1y_iii, v2x_iii, v2y_iii = run_Verlet(r1_iii, r2_iii) 
plt.figure(figsize=(8,6))
plt.plot(r1x_iii, r1y_iii, '.', label="particle 1")
plt.plot(r2x_iii, r2y_iii, '.', label="particle 2" )
plt.legend(fontsize=14)
plt.xlabel("x(t)", fontsize=14)
plt.ylabel("y(t)", fontsize=14)
plt.title("Trajectories of particles for the initial condition (iii)", fontsize=16)
plt.savefig("2c.pdf")


h = 0.01 # dt
num_steps = 100
t = np.linspace(0, h*(num_steps-1),num_steps)

plt.figure(figsize=(8,6))
plt.plot(t, r1x_i, '.', label="Particle1")
plt.plot(t, r2x_i, '.', label="Particle2")
plt.xlabel("time (s)", fontsize=14)
plt.ylabel("x-position x(t)", fontsize=14)
plt.legend(fontsize=14)
plt.title("X position vs. time", fontsize=16)
plt.savefig("2d.pdf")



    
