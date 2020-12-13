# Author: Aslesha Pokhrel

# Q2
import numpy as np
import matplotlib.pyplot as plt
from pylab import clf, plot, pause, draw

L = 1.0
g = 9.81
n_b=0.0
H = 0.01
J = 51
dx = 0.02
x = np.arange(0, L+dx, dx)
dt = 0.01
A = 0.002
mu = 0.5
sigma = 0.05
tmax = 5
eps = 1e-9

u = np.zeros(x.shape, dtype=float)
gauss =  A*np.exp(- (x - mu)**2 / (sigma**2))
eta = H + gauss - np.mean(gauss)

u_new = np.zeros(x.shape, dtype=float)
eta_new = np.zeros(x.shape, dtype=float)


t = 0.0
while t < tmax:
# for t in t_arr:
    
    for j in range(0, J):
        if j == 0:
            u_new[0] = 0.0
            eta_new[j] = eta[j] - dt / dx * (eta[1] * u[1] - eta[0]*u[0])
        elif j == J-1:
            u_new[j] = 0.0
            eta_new[j] = eta[j] - dt / dx * (eta[j] * u[j] - eta[j-1] *u[j-1])
        else:
            Fj_plus1 = 0.5 * (u[j+1])**2  + g * eta[j+1]
            Fj_minus1 = 0.5 * (u[j-1])**2 + g * eta[j-1]
            u_new[j] = u[j] - dt / (2*dx) * (Fj_plus1 - Fj_minus1)         
            eta_new[j] = eta[j] - dt / (2*dx) * (eta[j+1] * u[j+1] - eta[j-1]*u[j-1])
            
    #         if j == 10:
    #             print(u_new[j] )
    #             break
    # break
            
    u = np.copy(u_new)
    eta = np.copy(eta_new)
    
    for ti in [0,1,4]:      
        if np.abs(t-ti)<eps:
            
            plt.plot(x, eta, label=f"$t =${round(t, 2)}s")
    
    # clf() # clear the plot
    # plot(x, eta) # plot the current sin curve
    # draw()
    # pause(0.001) #pause to allow a smooth animation
    
    t += dt

plt.legend()






