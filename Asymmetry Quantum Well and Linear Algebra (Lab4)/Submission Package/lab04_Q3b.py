import numpy as np
import matplotlib.pyplot as plt

# set the value of c
c = 2

# set the initial value of x
x0 = 1

# desired accuracy
accuracy = 1e-6

# initialize difference in x to 1
dx = 1

# define the function of interest
def f(x):
    return 1 - np.exp(-c*x)

# list to store subsequent x0
x0_list = [x0]

# set the iteration count to 0
iteration_count = 0

while dx > accuracy:
    # increment the iteratioon count
    iteration_count += 1
    # append the new value to the list
    x0_list.append(f(x0_list[-1]))
    # set the dx to the difference of last two x0s
    dx = np.abs(x0_list[-1] - x0_list[-2])
    
# print the number of iterations and final solution
print("Using the relaxation method:")
print("\tIt takes {} iterations to converge to a " 
        "solution accurate to 1e-6".format(iteration_count))
print("\tIt converges to x = {}".format(x0_list[-1]))


# Part(c) overrelaxation  method to find the root

# define new function for overrelaxation method
def g(omega, x):
    # use f(x) from above
    return (1 + omega) * f(x) - omega * x


# set the value of omega
omega = np.arange(0, 1, 0.1)

# store the number of iterations for each omega
iter_count_list = []

# store the final value of x for each omega
final_x_list = []

# try different omega between 0 and 1
for w in omega:
    # set the iteration count to 0
    iteration_count_overrelaxation = 0
    
    # list to store subsequent x0
    x0_overrelaxation_list = [x0]
    
    # initialize difference in x to 1
    dx = 1
    
    while dx > accuracy:
        # increment the iteratioon count
        iteration_count_overrelaxation += 1
        # append the new value to the list
        x0_overrelaxation_list.append(g(w, x0_overrelaxation_list[-1]))
        # set the dx to the difference of last two x0s
        dx = np.abs(x0_overrelaxation_list[-1] - x0_overrelaxation_list[-2])
    
    # store the final x value
    final_x_list.append(x0_overrelaxation_list[-1])
    # store the number of iterations
    iter_count_list.append(iteration_count_overrelaxation)


# plot to see the optimal omega
plt.figure()
plt.plot(omega, iter_count_list)
plt.xlabel("omega value")
plt.ylabel("num iterations")
plt.title("Number of iterations vs. omega")


# print the number of iterations
print("\nFor overrelaxation method, omega = 0.5 gives the lowest number of iterations to"
      " converge to a solution accurate to 1e-6.")
print("\t In this case, the number of required iteration is {}".format(min(iter_count_list)))    

# find the corresponding index of x 
ind = iter_count_list.index(min(iter_count_list))

# print the answer
print("\t And, the value of x is {}".format(final_x_list[ind]))

    