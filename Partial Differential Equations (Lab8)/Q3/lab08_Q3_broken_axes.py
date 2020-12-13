# Author : Aslesha Pokhrel

"""
Disclaimer: Run "pip install brokenaxes" in the command prompt before 
running this file. If you are not willing to do so, run lab08_Q3.py only.
"""

import numpy as np
import matplotlib.pyplot as plt
# for q3b plot
from brokenaxes import brokenaxes
from pylab import clf, plot, pause, draw


def F_u(u, eta, g=9.81):
    """
    Return first element of F(u, eta)
    """
    return 0.5 * (u)**2  + g * eta


def F_eta(u, eta, eta_b=0):
    """
    Return second element of F(u, eta).
    """
    return (eta - eta_b) * u


def main(q3a=True, animation=False):
    
    # setting the given conditions
    L = 1.0
    H = 0.01
    
    if q3a:
        J = 50 # q3a
        dx = 0.02 # q3a
        dt = 0.01 # q3a
        A = 0.002 # q3a
        sigma = 0.05 # q3a
        mu = 0.5 # q3a
    else:
        J=150 # q3b
        dx = 1/150 # q3b
        dt = 0.001 # q3b
        A = 2 * 10**(-4) # q3b
        sigma = 0.1 #q3b
    
    x = np.arange(0, L, dx)
    # end time 
    tmax = 5
    # start time
    t = 0.0
    # eps to check difference in time 
    eps = 1e-9
        
    # set initial conditions
    u = np.zeros(x.shape, dtype=float)
    
    # initialization
    if q3a:
        gauss =  A*np.exp(- (x - mu)**2 / (sigma**2)) # q3a
        eta = H + gauss - np.mean(gauss) # q3a
        eta_b_arr = np.zeros(x.shape, dtype=float) # q3a 
    else:
        gauss =  A*np.exp(-x**2 / (sigma**2)) # q3b
        eta = H + gauss - np.mean(gauss) # q3b
        eta_bs = H - (4*10**(-4)) # q3b
        alpha = (8*np.pi) # q3b
        x0 = 0.5 # q3b
        eta_b_arr = (1 + np.tanh(alpha * (x - x0))) * eta_bs/2  # q3b
    
    # more arrays to hold the results
    u_new = np.zeros(x.shape, dtype=float)
    eta_new = np.zeros(x.shape, dtype=float)
    
    u_half = np.zeros(x.shape, dtype=float)
    eta_half = np.zeros(x.shape, dtype=float)
    
    # to plot figure
    plt.figure(figsize=(8,6))
    if not q3a:
        # plot broken axis for q3b
        bax = brokenaxes(ylims=((-0.0001, .0097), (.0099, .0102)), hspace=.03, 
                         height_ratios=[2,1])
    # constant that comes up frequently
    const = dt / dx 
    while t < tmax:
        for j in range(0, J):
            if j == 0:
                # boundary condition
                u_new[0] = 0.0 
                # handle boundry points with forward difference
                eta_new[j] = eta[j] - const * (F_eta(eta[1], u[1], eta_b_arr[1])
                                                    - F_eta(eta[0], u[0], eta_b_arr[0]))
                # calculate u and eta at half-step eq.12
                u_half[j] = 0.5 * (u[j+1] + u[j]) - const/2 * (F_u(u[j+1], eta[j+1])
                                                              - F_u(u[j], eta[j]))
                eta_half[j] = 0.5 * (eta[j+1] + eta[j]) - const/2 * (F_eta(u[j+1], eta[j+1], eta_b_arr[j+1]) 
                                                                    - F_eta(u[j], eta[j], eta_b_arr[j]))
            elif j == J-1:
                # boundary condition
                u_new[j] = 0.0
                # handle boundry points with forward difference
                eta_new[j] = eta[j] - const * (F_eta(eta[j], u[j], eta_b_arr[j])
                                                    - F_eta(eta[j-1], u[j-1], eta_b_arr[j-1]))
            else:   
                # calculate u and eta at half-step eq.12
                u_half[j] = 0.5 * (u[j+1] + u[j]) - const/2 * (F_u(u[j+1], eta[j+1])
                                                              - F_u(u[j], eta[j]))
                eta_half[j] = 0.5 * (eta[j+1] + eta[j]) - const/2 * (F_eta(u[j+1], eta[j+1], eta_b_arr[j+1])
                                                                    - F_eta(u[j], eta[j], eta_b_arr[j]))
                # update u and eta eq.13
                u_new[j] = u[j] - const * (F_u(u_half[j], eta_half[j]) 
                                              - F_u(u_half[j-1], eta_half[j-1]))        
                eta_new[j] = eta[j] - const * (F_eta(u_half[j], eta_half[j], (eta_b_arr[j] + eta_b_arr[j-1])/2 ) 
                                                  - F_eta(u_half[j-1], eta_half[j-1], (eta_b_arr[j] + eta_b_arr[j-1])/2 ))  
        # update u and eta       
        u = np.copy(u_new)
        eta = np.copy(eta_new)
        
        if animation:
            clf() # clear the plot
            plot(x, eta) # plot the current sin curve
            draw()
            pause(0.001) #pause to allow a smooth animation
        else:
            # plot at certain t for 3a and 3b 
            if q3a:
                twanted = [0,1,4]
            else:
                twanted = [0,1,2,4]
            for ti in twanted:      
                if np.abs(t-ti)<eps:
                    if q3a:
                        # use matplotlib for 3a
                        plt.plot(x, eta, label=f"$t =${round(t, 2)}s")
                    else:
                        # use brokenaxes for 3b
                        bax.plot(x, eta, label=f"$t =${round(t, 2)}s")
        # increment time
        t += dt
    
    if not animation:
        # set the plotting stuff
        if q3a:
            plt.legend(fontsize=12)
            plt.xlabel("Position $x$", fontsize=14)
            plt.ylabel("Altitiude of free surface $\eta(x,t)$", fontsize=14)
            plt.title("Simulation of the shallow water system with Two-Step Lax-Wendroff scheme")
            plt.xticks(fontsize=12) # change size of the x-axis ticks
            plt.yticks(fontsize=12) # change size of the y-axis ticks
            plt.tight_layout()
            plt.savefig("q3a.pdf")
            plt.show()
        else:
            # plot bottom topography for q3b
            bax.plot(x, eta_b_arr, label="bottom topography")
            bax.set_xlabel("Position $x$", fontsize=14)
            bax.set_ylabel("Altitiude ($\eta$)", fontsize=14)
            bax.legend()
            plt.title("Tsunami simulation with variable bottom topography")
            plt.savefig("q3b.pdf")
            plt.show()

# q3a
main()
# q3b
main(False, True)
