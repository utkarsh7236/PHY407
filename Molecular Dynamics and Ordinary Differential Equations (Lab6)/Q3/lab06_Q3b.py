import numpy as np
import matplotlib.pyplot as plt
from lab06_Q3_functions import distance
from lab06_Q3ac import run_verlet
import seaborn as sns

# for formatting the graph
sns.set()


def KE(vx, vy, m=1):
    """
    Return the kintetic energy of a particle given its x and y velocity.
    """
    return 1/2 * m * (vx**2 + vy**2)


def PE(r, eps = 1, sigma = 1):
    """
    Return the potential energy per particle.
    """
    # Lennard - Jones Potential
    V = 4 * eps * ((sigma * r**(-1))**12 - (sigma * r**(-1))**6)
    # Potential energy per particle is given by V(r)/2
    return V/2  


# Q3 b
def get_energies(v, r, num_particles, num_steps):
    """
    Given the velocity and position of each particle for num_steps, return the
    potential, kinetic and total energy.
    """
    # initialize arrays to store energies at each time step
    PE_arr = np.zeros(num_steps) # potential energy
    KE_arr = np.zeros(num_steps) # kinetic energy
    TE_arr = np.zeros(num_steps) # total energy
    
    for t in range(num_steps):
        for i in range(num_particles):
            # need a second for loop for PE to calculate interaction of i with 
            # all j particles
            for j in range(num_particles):
                if i != j:
                    PE_arr[t] += PE(distance(r[i, t, 0], r[i, t, 1], r[j, t, 0], 
                                             r[j, t, 1])) # calculate and accumulate PE for each particle
        
            KE_arr[t] += KE(v[i, t, 0], v[i, t, 1]) # accumulate KE of ith particle
            
        TE_arr = KE_arr + PE_arr
        
    return PE_arr, KE_arr, TE_arr


def plot_result(is_periodic=False):
    # run the program
    num_particles = 16
    num_steps = 1000
    h = 0.01
    # get velocity and postion arrays
    v, r = run_verlet(num_particles, num_steps, h, is_periodic)
    
    # calculate energies
    PE_arr, KE_arr, TE_arr = get_energies(v, r, num_particles, num_steps)
    
    ############# Plot the energies #############
    # initiaize time array
    t_arr = np.linspace(0,h*(num_steps-1),num_steps) 

    plt.figure(figsize=(9,8))
    plt.plot(t_arr, TE_arr, '.', label='Total')
    plt.plot(t_arr, KE_arr, '.', label='Kinetic')
    plt.plot(t_arr, PE_arr, '.', label='Potential')
    plt.legend(fontsize=14)
    plt.xlabel("time (s)", fontsize=16)
    plt.ylabel("Energy", fontsize=16)
    title = "Potential, Kinetic and Total energy of the particles vs time"
    plt.title(title, fontsize=16)
    plt.tight_layout()
    plt.savefig("3bb.pdf")
    
    # plot total energy within -1, +1
    plt.figure(figsize=(10,8))
    plt.plot(t_arr, TE_arr, '.', label='Total')
    plt.plot(t_arr, np.ones(num_steps)* np.mean(TE_arr)*1.01, '.', label='mean(TE) - 1%')
    plt.plot(t_arr, np.ones(num_steps)* np.mean(TE_arr)*0.99, '.', label='mean(TE) + 1%')
    plt.legend(bbox_to_anchor=(0.90, 0.90), loc='upper right',fontsize=14)
    plt.xlabel("time (s)", fontsize=16)
    plt.ylabel("Total Energy", fontsize=16)
    title = "Total energy of the particles over time"
    plt.title(title, fontsize=16)
    plt.tight_layout()
    plt.savefig("3b.pdf")
    

if __name__ == "__main__":
    # plot the energies
    plot_result()

    
