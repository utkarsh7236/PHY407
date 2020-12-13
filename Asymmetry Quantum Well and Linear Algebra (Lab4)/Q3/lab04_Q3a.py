import numpy as np
import matplotlib.pyplot as plt

# set the value of c
c_list = np.arange(0.01, 3+0.01, 0.01)

# set the initial value of x
x0 = 1

# desired accuracy
accuracy = 1e-6

# define the function of interest
def f(c, x):
    return 1 - np.exp(-c*x)


# store the final x for each c
final_ans = []

# calculate x for all c in c_list
for c in c_list:
    # initialize difference in x to 1
    dx = 1
    
    # list to store subsequent x0
    x0_list = [x0]
    
    while dx > accuracy:
        # append the new value to the list
        x0_list.append(f(c, x0_list[-1]))
        # set the dx to the difference of last two x0s
        dx = np.abs(x0_list[-1] - x0_list[-2])
    
    # store final answer for each c in final_ans list
    final_ans.append(x0_list[-1])
    

plt.figure(figsize=(8,6))
plt.scatter(c_list, final_ans, s=7)
plt.xlabel("Constant (c)")
plt.ylabel("Root (x)")
plt.grid()
plt.title("Root of $x = 1 - e^{-cx}$ for various values of c")
plt.savefig("Q3a.pdf")

