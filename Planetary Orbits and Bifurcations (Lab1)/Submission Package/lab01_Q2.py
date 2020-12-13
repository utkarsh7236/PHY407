import numpy as np
import matplotlib.pyplot as plt


def calculate_population(initial_norm_pop, reproduction_rate, max_years):
    """
    initial_norm_pop : initial normalized population x0
    reproduction_rate : maximumreproduction rate r
    max_years : total number of years pmax 
    
    Return a list of values normalized population (xp), one for each year p.
    """
    population = [initial_norm_pop]
    for i in range(1, max_years):
        # logistic map eqn
        population.append(reproduction_rate * (1 - population[i-1]) * population[i-1])
    
    return population


def q2c():
    """
    Plot the evolution of population for different values of r between 2 and 4 
    with max_years = 50 and initial population = 0.1
    """
    x0 = 0.1
    max_years = 50
    years_arr = np.arange(max_years)
    
    rates = [2, 2.5, 3, 3.5, 4] 
    
    plt.figure(figsize=(10,8))
    
    for r in rates:
        population = calculate_population(x0, r, max_years)
        plt.plot(years_arr, population, label="rate=" + str(r))
        
    plt.grid()
    plt.xlabel("years")
    plt.ylabel("normalized population (x_p)")
    plt.title("Evolution of population over years")
    plt.legend()
    plt.savefig("2c.pdf")
    

def q2de(rates, fig_name):
    """
    Plot bifurcation diagram.
    """
    x0 = 0.1
    max_years = 2000
    
    plt.figure()
    
    for r in rates:
        population = calculate_population(x0, r, max_years)
        if r < 3:
            r_arr = [r] * 100
            population = population[-100 : ]
        else:
            r_arr = [r] * 1000
            population = population[-1000 : ]
            
        plt.plot(r_arr, population, 'k.', markersize=0.1)
        
    plt.xlabel("rate (r)")
    plt.ylabel("normalized population (x_p)")
    plt.title("Bifurcation Diagram")
    plt.savefig(fig_name)
    

def q2fg():
    """
    Plot the evolution of population for perturbed initial condition and it's 
    difference with the non-perturbed initial condition.
    """
    x0 = 0.2
    eps = 10**(-13)
    r = 3.745
    max_years = 150
    years_arr = np.arange(max_years)
    
    pop1 = calculate_population(x0, r, max_years)
    pop2 = calculate_population(x0 + eps, r, max_years)
    
    plt.figure(figsize=(8,6))
    plt.plot(years_arr, pop1, label="x0")
    plt.plot(years_arr, pop2, label="x0+eps")
    plt.xlabel("years")
    plt.ylabel("normalized population (x_p)")
    plt.title("Evolution of population over years for perturbed intial condition")
    plt.legend(loc="lower right")
    plt.savefig("q2f.pdf")
    
    
    # q2(g)
    curve_fit = np.polyfit(years_arr[:100], np.log(np.abs(np.array(pop2) - np.array(pop1)))[:100], 1)
    print(curve_fit)
    
    lam, a = curve_fit[0], curve_fit[1]
    
    plt.figure()
    plt.semilogy(years_arr, np.abs(np.array(pop2) - np.array(pop1)), label="absolute difference")
    plt.plot(years_arr[:100], np.e**a * np.e ** (lam * years_arr[:100]), label="curve fit")
    plt.xlabel("years")
    plt.ylabel("absolute difference in population")
    plt.title("Absolute difference in population for perturbed intial condition")
    plt.legend()
    plt.savefig("q2g.pdf")
    

if __name__ == "__main__":
    q2c()
    q2de(np.arange(2, 4, 0.005), "q2d.png")
    q2de(np.arange(3.738, 3.745, 10**(-5)), "q2e.png")  
    q2fg()
    
    


