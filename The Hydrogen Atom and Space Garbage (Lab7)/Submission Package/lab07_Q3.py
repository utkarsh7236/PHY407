import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.constants as pc
from scipy.integrate import simps

# for formatting the graph
sns.set()

# set up the constants. From handout
a = pc.physical_constants['Bohr radius'][0]
E0 = pc.physical_constants['Rydberg constant times hc in eV'][0]
e_m = pc.m_e # electron mass
hbar = pc.hbar
e = pc.e # elementary charge
eps0 = pc.epsilon_0

# set h and r_infty
h = 0.0001 * a # step size
r_infty = 20 * a # max value of r

############# Code adapted from squarewell.py online #############
# Potential function
def V(x):
    """
    Return the potential given by equation 3.
    """
    return - e**2 / (4 * np.pi * eps0 * x)

def f(r,x,E,l):
    """
    Return the second order ODE in r as a pair of coupled first order ODEs for 
    R and S. 
    """
    # set R and S to given values in r
    R = r[0]
    S = r[1]
    
    # update R and S written as two first order ODEs from equation 2
    fR = S
    fS = l*(l+1)*R/x**2 + (2*e_m*R*(V(x)-E))/hbar**2 - 2*S/x
    
    # return the updated value
    return np.array([fR,fS],float)

# Calculate the wavefunction for a particular energy
def solve(E, l):
    R = 0.0
    S = 1.0
    r = np.array([R,S],float)
    
    # to store all the values of the R function
    R_arr = []
    
    # use RK4 to solve the ODE
    for x in np.arange(h,r_infty,h):
        # add the R value to R_arr
        R_arr.append(r[0])
        
        # update using RK4 method
        k1 = h*f(r,x,E,l)
        k2 = h*f(r+0.5*k1,x+0.5*h,E,l)
        k3 = h*f(r+0.5*k2,x+0.5*h,E,l)
        k4 = h*f(r+k3,x+h,E,l)
        r += (k1+2*k2+2*k3+k4)/6

    return r[0], np.array(R_arr)

##############################################################

def find_E(n, l):
    """
    Main program to find the energy using the secant method
    """
    # bracket the energies
    E1 = -15*e/n**2
    E2 = -13*e/n**2
    R2, R_arr = solve(E1, l)
    
    # set the target
    target = e/1000
    # find energy by using the secant method
    while abs(E1-E2)>target:
        R1 = R2
        R2, R_arr = solve(E2, l)
        E1, E2 = E2,E2-R2*(E2-E1)/(R2-R1)
    
    # return energy and eigenfunction R
    return E2/e, R_arr


def normalize(R):
    """
    Normalize the eigenfunction R by using the integral(|R(r)|^2).
    """
    # get the integral of |R1|^2 using simpsons rule
    norm = simps(np.abs(R)**2)
    # return by dividing R by the sqrt of norm
    return R/np.sqrt(norm)   


def plot_eigenfunc(numeric_R, analytic_R, n, l, fig_name):
    """
    Plot the normalized numerical and analytic eigenfunctions.
    """
    # normalize numerical R eigenfunction
    n_R = normalize(numeric_R)
    # normalize analytic R 
    a_R = normalize(analytic_R)
    
    # plot the normalized eigenfuctions
    plt.figure(figsize=(8,6))
    plt.plot(r, n_R, label="numerical")
    plt.plot(r, a_R, label="analytic")
    plt.title("Normalized eigenfunction R(r) for n={}, l={}".format(n,l), 
              fontsize=14)
    plt.xlabel("radius r (m)", fontsize=14)
    plt.ylabel("normalized R(r)", fontsize=14)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_name + ".pdf")
    plt.show()


##### Part B (find energies) #### 
# first case
n1, l1 = 1, 0
Energy1, R1 = find_E(n1, l1)
print("The energy for n={} and l={} is {} eV.".format(n1, l1, Energy1))

# second case
n2, l2 = 2, 0
Energy2, R2 = find_E(n2, l2 )
print("The energy for n={} and l={} is {} eV.".format(n2, l2 , Energy2))

# third case
n3, l3 = 2, 1
Energy3, R3 = find_E(n3, l3)
print("The energy for n={} and l={} is {} eV.".format(n3, l3, Energy3))

print("|E3 -E2| = ", Energy3 - Energy2)


#### Part D ####
# initialize the values of r for plotting
r = np.arange(h,r_infty,h)

##### first case #####
# analytic solution from the given link
analytic_R1 = np.exp(-r/a) * 2/(a**(3/2))
# plot the normalized eigenfuctions
plot_eigenfunc(R1, analytic_R1, n1, l1, '3a')


##### second case #####
# analytic solution from the given link
analytic_R2 = np.exp(-r/(2*a)) / (2*np.sqrt(2) * a**(3/2)) * (2- r/a) 
# plot the normalized eigenfuctions
plot_eigenfunc(R2, analytic_R2, n2, l2, '3b')


# third case
# analytic solution from the given link
analytic_R3 = np.exp(-r/(2*a)) / (2*np.sqrt(6) * a**(3/2)) * (r/a) 
# plot the normalized eigenfuctions
plot_eigenfunc(R3, analytic_R3, n3, l3, '3c')


