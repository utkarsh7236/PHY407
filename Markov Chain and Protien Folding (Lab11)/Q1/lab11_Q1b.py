from numpy import cos, pi, exp, sqrt
import numpy as np
from numpy.random import normal, seed
from random import random
import matplotlib.pyplot as plt
from matplotlib import rc

font = {'family': 'DejaVu Sans', 'size': 14}  # adjust fonts
rc('font', **font)

seed(2)

def f(x, y):
    """
    Given function 1b(i) to find the minimum.
    """
    return x**2 - cos(4*x*pi) + (y-1)**2


def g(x,y):
    """
    Given function 1b(ii) to find the minimum.
    """
    return cos(x) + cos(sqrt(2) * x) + cos(sqrt(3) * x) + (y-1)**2 


def run_b(func, q_name):
    # set the variables
    Tmax = 1.0
    Tmin = 1e-5
    tau = 1e4
    
    # Main loop
    t = 0
    T = Tmax
    # starting point
    x,y = 2,2
    # std and mean for gaussian distribution
    mean, std = 0, 1
    
    # to store the results
    time = []
    func_arr = []
    x_y_arr = []
    
    # compute the results intially
    time.append(t)
    x_y_arr.append([x,y])
    func_arr.append(func(x,y))
    
    while T>Tmin:
    
        # Cooling
        t += 1
        T = Tmax*exp(-t/tau)
        
        # sample from a gaussian dstribution
        dx, dy = normal(mean, std, 2) 
        # Monte Carlo moves
        new_x, new_y = x+dx, y+dy
        
        # find the change in energy
        new_func = func(new_x, new_y)
        delta_func = new_func - func_arr[-1]
        
        if q_name == 'bi':
            if random() < exp(-delta_func/T):
                # make move
                x, y = new_x, new_y
                # store the energy, time and moves
                time.append(t)
                x_y_arr.append([x,y])
                func_arr.append(func(x,y))
        else:
            if 0 < new_x < 50 and -20 < new_y < 20 and random() < exp(-delta_func/T):
                # make move
                x, y = new_x, new_y
                # store the energy, time and moves
                time.append(t)
                x_y_arr.append([x,y])
                func_arr.append(func(x,y))
    
    print("The value of (x, y) for Qb({}) is ({:.3f}, {:.3f}).".format(q_name[1:],x,y))
    
    x_y_arr = np.array(x_y_arr)
    
    # plot x as a function of time
    plt.figure()
    plt.scatter(time, x_y_arr[:, 0])
    plt.xlabel("Time (t)")
    plt.ylabel("x")
    plt.title("Qb({}) x as a function of time".format(q_name[1:]))
    plt.grid()
    plt.tight_layout()
    plt.savefig(q_name + '_x.pdf')
    
    # plot y as a function of time
    plt.figure()
    plt.scatter(time, x_y_arr[:, 1], color='#f97306')
    plt.xlabel("Time (t)")
    plt.ylabel("y")
    plt.title("Qb({}) y as a function of time".format(q_name[1:]))
    plt.grid()
    plt.tight_layout()
    plt.savefig(q_name + '_y.pdf')
    
    
# run b(i)
run_b(f, 'bi')

# run b(ii)
run_b(g, 'bii')   
    
