"""
Strter code for protein folding
Author: Nicolas Grisuard, based on a script by Paul Kushner
"""

from random import random, randrange
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def calc_energy(monomer_coords, monomer_array):
    """ Compute energy of tertiary structure of protein """
    energy = 0.0

    # compute energy due to all adjacencies (incl. directly bonded monomers)
    for i in range(N):
        for nghbr in [[-1, 0], [1, 0], [0, -1], [0, 1]]:  # 4 neighbours
            nghbr_monomer = monomer_array[monomer_coords[i, 0] + nghbr[0],
                                          monomer_coords[i, 1]+nghbr[1]]

            if nghbr_monomer == 1:  # check neighbour is not empty
                energy += eps

    # divide by 2 to correct for double-counting
    energy = .5*energy

    # correct energy to not count directly bonded monomer neighbours
    energy -= (N-1)*eps

    return energy


def dist(position1, position2):
    """ Compute distance """
    return ((position1[0]-position2[0])**2+(position1[1]-position2[1])**2)**.5


font = {'family': 'DejaVu Sans', 'size': 14}  # adjust fonts
rc('font', **font)
dpi = 150

eps = -5.0  # interaction energy
N = 30  # length of protein
T_max = 10  # starting temperature for Monte Carlo
T_min = 0.5  # ending temperature for Monte Carlo
n = 5*int(1e5)  # number of Monte Carlo steps
dt = 0.5 # change in temp

# to store the results
temp_arr = []
energy_mean_arr = []
energy_std_arr = []

while T_max >= T_min: # run simulation starting from T=10 until T=0.5

    energy_array = np.zeros(n)  # initialize array to hold energy
    
    # initialize arrays to store protein information
    # 1st column is x coordinates, 2nd column is y coordinates, of all N monomers
    monomer_coords = np.zeros((N, 2), dtype='int')
    
    # initialize position of polymer as horizontal line in middle of domain
    monomer_coords[:, 0] = range(N//2, 3*N//2)
    monomer_coords[:, 1] = N
    
    # 2D array representing lattice,
    # equal to 0 when a lattice point is empty,
    # and equal to 1 when there is a monomer at the lattice point
    monomer_array = np.zeros((2*N+1, 2*N+1), dtype='int')
    
    # fill lattice array
    for i in range(N):
        monomer_array[monomer_coords[i, 0], monomer_coords[i, 1]] = 1
    
    # calculate energy of initial protein structure
    energy = calc_energy(monomer_coords, monomer_array)
    
    # do Monte Carlo procedure to find optimal protein structure
    for j in range(n):
        energy_array[j] = energy
    
        # move protein back to centre of array
        shift_x = int(np.mean(monomer_coords[:, 0])-N)
        shift_y = int(np.mean(monomer_coords[:, 1])-N)
        monomer_coords[:, 0] -= shift_x
        monomer_coords[:, 1] -= shift_y
        monomer_array = np.roll(monomer_array, -shift_x, axis=0)
        monomer_array = np.roll(monomer_array, -shift_y, axis=1)
    
        # pick random monomer
        i = randrange(N)
        cur_monomer_pos = monomer_coords[i, :]
    
        # pick random diagonal neighbour for monomer
        direction = randrange(4)
    
        if direction == 0:
            neighbour = np.array([-1, -1])  # left/down
        elif direction == 1:
            neighbour = np.array([-1, 1])  # left/up
        elif direction == 2:
            neighbour = np.array([1, 1])  # right/up
        elif direction == 3:
            neighbour = np.array([1, -1])  # right/down
    
        new_monomer_pos = cur_monomer_pos + neighbour
    
        # check if neighbour lattice point is empty
        if monomer_array[new_monomer_pos[0], new_monomer_pos[1]] == 0:
            # check if it is possible to move monomer to new position without
            # stretching chain
            distance_okay = False
            if i == 0:
                if dist(new_monomer_pos, monomer_coords[i+1, :]) < 1.1:
                    distance_okay = True
            elif i == N-1:
                if dist(new_monomer_pos, monomer_coords[i-1, :]) < 1.1:
                    distance_okay = True
            else:
                if dist(new_monomer_pos, monomer_coords[i-1, :]) < 1.1 \
                   and dist(new_monomer_pos, monomer_coords[i+1, :]) < 1.1:
                    distance_okay = True
    
            if distance_okay:
                # calculate new energy
                new_monomer_coords = np.copy(monomer_coords)
                new_monomer_coords[i, :] = new_monomer_pos
    
                new_monomer_array = np.copy(monomer_array)
                new_monomer_array[cur_monomer_pos[0], cur_monomer_pos[1]] = 0
                new_monomer_array[new_monomer_pos[0], new_monomer_pos[1]] = 1
    
                new_energy = calc_energy(new_monomer_coords, new_monomer_array)
    
                if random() < np.exp(-(new_energy-energy)/T_max):
                    # make switch
                    energy = new_energy
                    monomer_coords = np.copy(new_monomer_coords)
                    monomer_array = np.copy(new_monomer_array)
    
    # store the results
    temp_arr.append(T_max)
    energy_mean_arr.append(np.mean(energy_array))
    energy_std_arr.append(np.std(energy_array))
    # step by dt for next round                
    T_max -= dt

# plot the mean energy for each temp with errorbar
plt.figure(figsize=(8,6))
plt.errorbar(temp_arr, energy_mean_arr, yerr=energy_std_arr, fmt="-o", capsize=3)
plt.xlabel("Temperature (T)")
plt.ylabel("Energy")
plt.title("Energy vs. Temperature plot with errorbar")
plt.savefig("q3d.pdf")
                
