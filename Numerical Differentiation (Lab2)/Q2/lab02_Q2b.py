import numpy as np
import matplotlib.pyplot as plt
from lab02_Q2_functions import simpson_integral
from scipy.special import jv


def f(m, x, theta):
    """
    Helper function to calculate function inside the integral in Bessel function
    i.e. cos(m*theta - x*sin(theta)).
    m, x : Values to evaluate Bessel function at.
    """
    return np.cos(m * theta - x * np.sin(theta))
    

def J(m, x):
    """
    Evaluate Bessel function using Simpson's rule.
    m, x : Values to evaluate Bessel function at.
    """
    # Number of slices to use in Simpson's integral
    N = 1000
    # use Simpson integral to integrate helper function f from 0 to pi. Using
    # lambda function to pass in static argumnents to f.
    integral_result = simpson_integral(N, lambda theta : f(m, x, theta), 0, np.pi)
    
    # Multiply by 1/pi and return Bessel function evaluation
    return integral_result / np.pi


def q2ba():
    """
    Program for running q2(b) part (a) or exercise 5.4 (a).
    """
    
    # Nested list to store the final results as [[J(0,0), ... J(0,20)], 
    # [J(1,0), ... J(1,20)], [J(2,0), ... J(2,20)]] produced by J(m,x)
    bessel_arr = [[], [], []] 
    
    # Nested list to store the final results in the same format as bassel_arr
    # produced by scipy.special.jv
    jv_arr = [[], [], []]
    
    for m in range(3): # for J_0, J_1, J_2
        for x in range(21): # for x = 0, 1, ..., 20
            # evaluate Bessel function using function J code above
            bessel_arr[m].append(J(m,x))
            
            # evlauate function using scipy.special.jv function
            jv_arr[m].append(jv(m,x))
    
    
    # x array from 0 to 20 to be used for the plot
    x_arr = np.arange(21) 
    
    # Plot Bessel function for m = 0,1,2 and x = 0,1,...,20 in the same plot
    # evaluated using J(m,x) function coded above
    plt.figure()
    plt.plot(x_arr, bessel_arr[0], label="J_0") # Plot J_0 against x
    plt.plot(x_arr, bessel_arr[1], label="J_1") # Plot J_1 against x
    plt.plot(x_arr, bessel_arr[2], label="J_2") # Plot J_2 against x
    plt.xlabel("x")
    plt.ylabel("J_m(x)")
    plt.title("Bessel functions as a function of x using Simpson's rule")
    plt.legend()
    plt.savefig("q2b_simp.png")
    
    # Plot Bessel function for m = 0,1,2 and x = 0,1,...,20 in the same plot
    # evaluated using scipy.special.jv function
    plt.figure()
    plt.plot(x_arr, jv_arr[0], label="J_0") # Plot J_0 against x
    plt.plot(x_arr, jv_arr[1], label="J_1") # Plot J_1 against x
    plt.plot(x_arr, jv_arr[2], label="J_2") # Plot J_2 against x
    plt.xlabel("x")
    plt.ylabel("scipy.special.jv(x)")
    plt.title("Bessel functions as a function of x using scipy.special.jv")
    plt.legend()
    plt.savefig("q2b_jv.png")    
    

def diffraction_intensity(wave_length, r):
    """
    Return the intensity of light in the diffraction pattern for telescope.
    wave_length : wave length of light
    r          : the distance in the focal plane from the center of the 
                  diffraction pattern
    """
    # k = (2*pi)/lambda
    k = (2 * np.pi) / wave_length
    
    # return the intensity given by (J_1(kr)/(kr))^2
    return (J(1, k*r) / (k*r)) ** 2


def q2bb():
    """
    Program for q2(b) part (b) or exercise 5.4 (b).
    """
    wave_length = 0.5 # wavelength of 500 nm converted to micro-m
    # Create a grid from 0 to 1 micro-m
    x, y = np.mgrid[-1:1:200j, -1:1:200j]
    # The radius is given by sqrt(x^2 + y^2)
    r = np.sqrt(x**2 + y**2)
    # create a diffraction array
    diffraction_arr = diffraction_intensity(wave_length, r)
    
    # create the image using "hot" scheme
    plt.figure()
    plt.imshow(diffraction_arr, cmap="hot", vmax=0.01, extent=(-1,1,-1,1))
    plt.savefig("q2bb.pdf")
    

if __name__ == "__main__":
    q2ba()
    q2bb()
    
