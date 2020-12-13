""" Pseudocode for Trapezoidal rule: """
# Get the function (f) to evaluate
# Divide the interval [a, b] into N slices
# The width of the slice is given by h = (b-a)/N
# Implement eq. (5.3) pg. 142 from the textbook 
#   Accumulate the starting and ending bits in a sum as: s = 0.5*f(a) + 0.5*(b)
#   Evaluate and add the interior bits to the sum as:
#       for k from 1 to N-1:
#           add f(a+ k*h) to s
#   Multiply s by h and return
    

def trapezoidal_integral(N, f, a, b):
    """
    Calculate integral of one-variable function f from a to b with N slices 
    using Trapezoidal rule. Adapted from textbook and lecture notes.
    N : number of slices
    f : function with one argument
    a : starting point
    b : ending point
    """
    # width of the slice
    h = (b-a)/N
    # the starting and ending bits
    s = 0.5*f(a) + 0.5*f(b)
    # Evaluate and add the interior points excluding the start and end points
    for k in range(1, N): 
        s += f(a + k * h) # adding the interior bits
    
    # multiply the sum by h and return
    return h * s 


""" Pseudocode for Simpson's rule: """
# Get the function (f) to evaluate
# Divide the interval [a, b] into N slices. The number of has to be even
# The width of the slice is given by h = (b-a)/N
# Implement eq. (5.9) pg. 146 from the textbook
#   Accumulate the starting and ending bits in a sum as: s = f(a) + f(b)
#   Evaluate and add the interior bits at the odd points multiplied by 4 to the sum:  
#       for odd k from 1 to N-1:
#           add 4 * f(a+k*h) to the sum
#   Evaluate and add the interior bits at the even points multiplied by 2 to the sum:  
#       for even k from 1 to N-1:
#           add 2 * f(a+k*h) to the sum    
#   Multiply s by (1/3)*h and return 


def simpson_integral(N, f, a, b):
    """
    Calculate integral of one-variable function f from a to b with N slices 
    using Simpson's rule.
    N : number of slices
    f : function with one argument
    a : starting point
    b : ending point
    """
    # width of the slice
    h = (b-a)/N
    # the starting and ending bits
    s = f(a) + f(b)
    # adding the interior bits for odd k
    for k in range(1, N, 2):
        s += 4 * f(a + k * h) # multiply odd k part by 4
        
    # adding the interior bits for even k
    for k in range(2, N, 2):
        s += 2 * f(a + k * h) # multiply even k part by 2
    
    # multiply the sum by 1/3*h and return
    return 1 / 3 * h * s