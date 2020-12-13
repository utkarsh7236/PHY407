import numpy as np
import matplotlib.pyplot as plt
from lab06_Q3_functions import distance, acceleration
import seaborn as sns

# for formatting the graph
sns.set()


def total_acceleration(particle_i_index, r, t, num_particles, is_periodic=False):
    """
    particle_i_index : Index of the particle i in interest
    r : 3 dimensional array of particle positions
    t : time to identify proper position
    Return the total acceleration on particle i at time t.
    """
    # Given length
    Lx = 4.0
    Ly = 4.0
    # initialize x and y comp of the acceleration on each particle i
    acc_x = 0 
    acc_y = 0
    
    if is_periodic:
        # array to store position of replicated particles in 8 domain
        replicated_particles = np.zeros([8, num_particles, 2]) 
        # replicate the position of the particle in the domain 8 times
        # upper left corner
        replicated_particles[0, :, 0] = r[:, t, 0] - Lx
        replicated_particles[0, :, 1] = r[:, t, 1] + Ly
        # top middle
        replicated_particles[1, :, 0] = r[:, t, 0] 
        replicated_particles[1, :, 1] = r[:, t, 1] + Ly
        # upper right
        replicated_particles[2, :, 0] = r[:, t, 0] + Lx
        replicated_particles[2, :, 1] = r[:, t, 1] + Ly
        # left
        replicated_particles[3, :, 0] = r[:, t, 0] - Lx
        replicated_particles[3, :, 1] = r[:, t, 1] 
        # right
        replicated_particles[4, :, 0] = r[:, t, 0] + Lx
        replicated_particles[4, :, 1] = r[:, t, 1] 
        # bottom left
        replicated_particles[5, :, 0] = r[:, t, 0] - Lx
        replicated_particles[5, :, 1] = r[:, t, 1] - Ly
        # bottom middle
        replicated_particles[6, :, 0] = r[:, t, 0]
        replicated_particles[6, :, 1] = r[:, t, 1] - Ly
        # bottom right
        replicated_particles[7, :, 0] = r[:, t, 0] + Lx
        replicated_particles[7, :, 1] = r[:, t, 1] - Ly
                
    # accumulate the acceleration due to all the interaction for particle i
    for j in range(num_particles):
        if particle_i_index != j: # only accumulate the acceleration if the particle is not i itself
            xi =  r[particle_i_index, t, 0]
            yi =  r[particle_i_index, t, 1]
            xj =  r[j, t, 0]
            yj =  r[j, t, 1]
            # multiply by  x-dispacement/distance(xi, yi, xj, yj) for x component
            acc_x += acceleration(xi, yi, xj, yj) * ((xj - xi) / 
                                  distance(xi, yi, xj, yj))
            # multiply by  y-dispacement/distance(xi, yi, xj, yj) for y component
            acc_y += acceleration(xi, yi, xj, yj) * ((yj - yi) / 
                                  distance(xi, yi, xj, yj))
            
            if is_periodic:
                # add up acceleration from the 8 replicated particles as well
                for m in range(8):
                    # position of particles from the replicated_particles
                    xj =  replicated_particles[m, j, 0]
                    yj =  replicated_particles[m, j, 1]

                    # for x-comp
                    acc_x += acceleration(xi, yi, xj, yj) * ((xj - xi) / 
                                          distance(xi, yi, xj, yj))
                    # for y-comp
                    acc_y += acceleration(xi, yi, xj, yj) * ((yj - yi) / 
                                          distance(xi, yi, xj, yj))
                    
    return acc_x, acc_y 


def run_verlet(num_particles, num_steps, h, is_periodic=False):
    """
    num_particles : number of particles in the simulation
    num_steps : total number of time steps T. Default to 1000.
    h : step size dt. Default to 0.01
    """
    
    # setting the initial positions of the particles to be evenly spaced in a 
    # square domain using the code snippet given in the handout
    Lx = 4.0
    Ly = 4.0
    dx = Lx/np.sqrt(num_particles)
    dy = Ly/np.sqrt(num_particles)
    x_grid = np.arange(dx/2, Lx, dx)
    y_grid = np.arange(dy/2, Ly, dy)
    xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
    x_initial = xx_grid.flatten()
    y_initial = yy_grid.flatten()
    
    # initialize arrays to hold the results
    v = np.zeros([num_particles, num_steps, 2]) # array for velocity. Initialize to 0
    r = np.zeros([num_particles, num_steps, 2]) # array for distance
    v_temp = np.zeros([num_particles, 2]) # to store v(t + h/2) in eq7 for each particle
    
    # set the initial position in r array
    r[:, 0, 0] = x_initial # update initial x pos
    r[:, 0, 1] = y_initial # update initial y pos
    
    # find v_temp in eq 7
    for i in range(num_particles):
        # use the helper function to get acceleration x-y components
        acc_x, acc_y =  total_acceleration(i, r, 0, num_particles, is_periodic)            
                
        # update x-component of v_temp for each particle
        v_temp[i, 0] = v[i, 0, 0] + 1/2 * h * acc_x
        # update y-component of v_temp for each particle
        v_temp[i, 1] = v[i, 0, 1] + 1/2 * h * acc_y
    
    # run Verlet method for num_steps time steps
    for t in range(1, num_steps):
        # calculate equation 8 of Verlet method
        r[:, t, 0] = r[:, t-1, 0] + h * v_temp[:, 0] # update x comp of distance r
        r[:, t, 1] = r[:, t-1, 1] + h * v_temp[:, 1] # update y comp of distance r
        
        # First change for the periodic boundary condition
        if is_periodic:
            # move particles inside the appropriate domain if they have 
            # exited the domain
            r[:, t, 0] = np.mod(r[:, t, 0], Lx)
            r[:, t, 1] = np.mod(r[:, t, 1], Ly)
            
                
        # calculate equation 9, 10 and 11
        for i in range(num_particles):
            # use the helper function to get acceleration x-y components
            acc_x, acc_y =  total_acceleration(i, r, t, num_particles, is_periodic)                        
            
            # calculate kx and ky of equation 9
            kx = h * acc_x
            ky = h * acc_y
            
            # calculate equation 10
            v[i, t, 0] = v_temp[i, 0] + 1/2 * kx # update x-comp of v
            v[i, t, 1] = v_temp[i, 1] + 1/2 * ky # update y-comp of v
            
            # calculate equation 11
            v_temp[i, 0] += kx # update x-comp of v(t + h/2)
            v_temp[i, 1] += ky # update y-comp of v(t + h/2)
    
    return v, r 


def plot_result(is_periodic=False):
    # run the program
    num_particles = 16
    num_steps = 1000
    h = 0.01
    # get velocity and postion arrays
    v, r = run_verlet(num_particles, num_steps, h, is_periodic)
    
    ## Plot the trajectory
    plt.figure(figsize=(10,8))
    for n in range(num_particles):
        plt.plot(r[n,:,0], r[n,:,1], '.', label="Particle {}".format(n+1))
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',  fontsize=14)
    plt.xlabel("x-position x(t)", fontsize=16)
    plt.ylabel("y-position y(t)", fontsize=16)
    title = "Trajectories of all 16 particles"
    if is_periodic:
        title += " with periodic boundary conditions" 
    plt.title(title, fontsize=16)
    plt.tight_layout()
    if not is_periodic:
        plt.savefig("3a.pdf")
    else:
       plt.savefig("3c.pdf") 
        

if __name__ == "__main__":
    # # without periodic boundary condition
    plot_result()
    
    # with periodic boundary condition
    plot_result(True)
    
    
