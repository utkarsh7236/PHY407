import numpy as np
import matplotlib.pyplot as plt
from lab03_Q1_gaussxw import gaussxwab


def f(T_a, t_h, u):
    """
    Function inside the integral in the probability equation.
    T_a : hourly temperature in degree celsius.
    t_h : snow surface age in hours.
    u : point to evaluate f at
    """
    # evaluate mu_bar at T_a, t_h
    result_u_bar = u_bar(T_a, t_h)
    # evaluate delta at T_a
    result_delta = delta(T_a)
    # return the function inside the integral in probability equation
    return np.exp(- (result_u_bar - u)**2 / (2 * result_delta**2))


def u_bar(T_a, t_h):
    """
    Return mean wind speed.
    T_a : hourly temperature in degree celsius.
    t_h : snow surface age in hours.
    """
    # u_bar formula from the question sheet
    return 11.2 + 0.365 * T_a + 0.00706 * T_a ** 2 + 0.9 * np.log(t_h)


def delta(T_a):
    """
    Return the standard deviation of the wind speed.
    T_a : hourly temperature in degree celsius.
    """
    # delta formula from the question sheet
    return 4.3 + 0.145 * T_a + 0.00196 * T_a**2


def probability(u_10, t_h, ta):
    """
    Return the probability of blowing snow in the Canadian Prairies.
    mu_10 : Avg hourly windspeed at a height of 10m
    T_a : hourly temperature in degree celsius.
    t_h : snow surface age in hours. 
    """
    # evaluate the integral
    # num sample points
    N = 100
    # starting point 
    a = 0
    # get the weights and sample points
    x, w = gaussxwab(N, a, u_10)
    # the loop of summing wp[k] * f(xp[k]) implemented using numpy to evaluate the integral
    integral_evaluation =  np.sum(np.array(w) * f(ta, t_h, np.array(x)))
    
    # evaluate delta at T_a
    result_delta = delta(ta)
    # return the probability by multiplying integral_evaluation by 1/(sqrt(2pi)* delta)
    return (1 / (np.sqrt(2 * np.pi) * result_delta)) * integral_evaluation
    

def main():
    # Avg hourly windspeed at a height of 10m
    u_10 = np.array([6, 8, 10])
    # snow surface age
    t_h = [24, 48, 72]
    # T_a = average hourly temperature in degree celsius 
    T_a = np.arange(-60, 40, 1)
    
    # colours for each value of u_10
    colours = ('r', 'g', 'b')
    # line format for each t_h
    lines = ('-', '--', ':')
    
    # Create a figure of figsize
    plt.figure(figsize=(7.5,6))
    plt.clf()
    # plot for all the combination of u_10 and t_h
    for (i_u_10, colour) in zip(u_10, colours): # Group plots for same t_h value by colour
        for (i_t_h, line) in zip(t_h, lines): # Group plots for same t_h value by line
            # line style in the plot
            plot_str = colour + line
            # calculate probability for all values of T_a
            p = [probability(i_u_10, i_t_h, ta) for ta in T_a]
            # Plot Probability as a function of T_a using plot_str style
            plt.plot(T_a, p, plot_str, 
                     label = "$t_h$={0}, {2}={1}".format(i_t_h, i_u_10, "$u_{10}$"))
            plt.legend()
            # label axes and title
            plt.xlabel("Temperature ($^{\circ}$C)")
            plt.ylabel("Probability")
            plt.title("Probability of blowing snow vs time for various $t_h$ and $u_{10}$")
            
    plt.savefig("q1b.pdf")
    
    
if __name__ == "__main__":
    main()
    
    