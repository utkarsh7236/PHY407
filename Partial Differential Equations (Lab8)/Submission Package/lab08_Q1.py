"""
Author: Aslesha Pokhrel
Code adaptaed from Newman's laplace.py
"""

from numpy import zeros, copy, max
import matplotlib.pyplot as plt
import numpy as np
from time import perf_counter


def calculate_potential(omega=0):
    """
    Calculate potential using the Gauss-Seidel method with overrelaxation. When 
    omega=0, overrelaxation is not used.
    """
    # Constants
    M = 100         # Grid squares on a side
    target = 1e-6   # Target accuracy
    
    # Create arrays to hold potential values
    phi = zeros([M+1,M+1],float)
    
    # Main loop
    delta = 1.0
    while delta>target:
        
        # make a copy of phi to calculate delta
        phi_copy = copy(phi)
        
        # Calculate new values of the potential
        for i in range(1, M): # boundaries never updated
            for j in range(1, M):
                if i >= 20 and i <= 80 and j == 20:
                    phi[i, j] = 1 # set phi to +1 for the first plate
                elif i >= 20 and i <= 80 and j == 80:
                    phi[i, j] = -1 # set phi to -1 for the second plate    
                else:
                    # update using Gauss-Seidel method
                    phi[i, j] = (1 + omega)/4 * (phi[i+1, j] + phi[i-1, j]
                    + phi[i, j+1] + phi[i, j-1]) - omega * phi[i, j]
        
        # Calculate maximum difference from old values
        delta = max(abs(phi-phi_copy))
    return phi


def contour_plot(phi, title, savename):
    """
    Make a contour plot of the potential phi.
    phi : 2D array of potential
    title : title of the plot
    savename : name to save the plot
    """
    # Make a plot
    plt.figure(figsize=(8, 6))
    plt.imshow(phi)
    plt.xlabel("x (cm)", fontsize=16)
    plt.ylabel("y (cm)", fontsize=16)
    plt.title(title, fontsize=14)
    # change the ticks so that label goes from -5 to 5 cm
    plt.yticks(np.linspace(0, 100, 11, dtype=int), 
               np.linspace(-5, 5, 11, dtype=int))
    plt.xticks(np.linspace(0, 100, 11, dtype=int), 
               np.linspace(-5, 5, 11, dtype=int))
    cbar = plt.colorbar()
    cbar.set_label('Potential ($V$)', fontsize=16)
    cbar.ax.tick_params(labelsize=14) # change size of the colorbar ticks
    plt.gray()
    plt.tight_layout()
    plt.xticks(fontsize=14) # change size of the x-axis ticks
    plt.yticks(fontsize=14) # change size of the y-axis ticks
    plt.savefig(savename)
    plt.show()
    

def stream_plot(phi, title, savename):
    """
    Make a stream plot of the electric field lines using the potential phi.
    phi : 2D array of potential
    title : title of the plot
    savename : name to save the plot
    """
    # make meshgrid 
    x = np.linspace(-5, 5, 101) # does not include zero
    y = np.linspace(-5, 5, 101)
    X, Y = np.meshgrid(x, y)
    
    # get electric field
    Ey, Ex = np.gradient(-phi, y, x) # careful about order
    
    # plot the stream plot
    fig = plt.figure(figsize=(8, 6))
    strm = plt.streamplot(X, Y, Ex, Ey, color=phi1, linewidth=2, cmap='autumn')
    cbar = fig.colorbar(strm.lines)
    cbar.set_label('Potential ($V$)', fontsize=16)
    cbar.ax.tick_params(labelsize=14) # change size of the colorbarticks
    plt.title(title, fontsize=14)
    plt.xlabel('$x$ (cm)', fontsize=16)
    plt.ylabel('$y$ (cm)', fontsize=16)
    plt.tight_layout()
    plt.xticks(fontsize=14) # change size of the x-axis ticks
    plt.yticks(fontsize=14) # change size of the y-axis ticks
    plt.savefig(savename)
    plt.show()
    

# without overrelaxation
s1 = perf_counter()
phi1 = calculate_potential()
t1 = perf_counter() - s1 # measure runtime
print("The time taken to calculate the electrostatic potential using "+
      "Gauss-Seidel without overrelaxation {:.2f}s".format(t1))
contour_plot(phi1, "Potential of Capacitor using Gauss-Seidel without overrelaxation",
              "q1a.pdf")
stream_plot(phi1, 'Electric field lines using Gauss-Seidel without overrelaxation',
            'q1a_stream.pdf')


# with relaxation w = 0.1
s2 = perf_counter() 
phi2 = calculate_potential(0.1)
t2 = perf_counter() - s2 # measure runtime
print("The time taken to calculate the electrostatic potential using"+
      " Gauss-Seidel with overrelaxation with omega=0.1 is {:.2f}s".format(t2))
contour_plot(phi2, "Potential of Capacitor using Gauss-Seidel with $\omega=0.1$",
              "q1b.pdf")
stream_plot(phi1, 'Electric field lines using Gauss-Seidel with $\omega=0.1$',
            'q1b_stream.pdf')


# with relaxation w = 0.5
s3 = perf_counter()
phi3 = calculate_potential(0.5)
t3 = perf_counter() - s3 # measure runtime
print("The time taken to calculate the electrostatic potential using"+
      " Gauss-Seidel with overrelaxation with omega=0.5 is {:.2f}s".format(t3))
contour_plot(phi3, "Potential of Capacitor using Gauss-Seidel with $\omega=0.5$",
              "q1c.pdf")
stream_plot(phi1, 'Electric field lines using Gauss-Seidel with $\omega=0.5$',
            'q1c_stream.pdf')



