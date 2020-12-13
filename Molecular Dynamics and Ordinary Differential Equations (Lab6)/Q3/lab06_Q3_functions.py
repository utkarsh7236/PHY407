import numpy as np


def distance(xi, yi, xj, yj):
    """
    Return the distance between i and j particles
    """
    return np.sqrt((xj - xi)**2 + (yj - yi)**2)


# acceleration found in part (a)
def acceleration(xi, yi, xj, yj):
    """
    Return the aceleration in r direction. 
    To get the x-component: multiply by x-dispacement/distance(xi, yi, xj, yj)
    To get the y-componet: multiply by y-displacement/distance(xi, yi, xj, yj)
    xi, yi : (x,y) coordinates of particle i 
    xj, yj : (x,y) coordinates of particle j
    """
    # the distance between the two particles
    r = distance(xi, yi, xj, yj)
    # use F = - d(V(r))/dr = ma where m = 1 to find a. Eps, sigma = 1 
    return  24 * ( r**(-7) - 2 * r**(-13))
