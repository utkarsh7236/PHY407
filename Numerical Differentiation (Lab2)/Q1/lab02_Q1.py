# Author: Utkarsh Mali
import numpy as np
import matplotlib.pyplot as plt


# Question 1a
def f(x):
    x_sq = x ** 2
    return np.exp(-x_sq)


# Creating h's list from 10^-16 to 10^0
h = []
for n in range(17):
    h.append(10 ** (-n))
h.reverse()
# print(f"h values are {h}")

# Setting initial value of x
x_0 = 0.5


# define forward difference
def forward_difference(x, h_step):
    x_ = x + h_step
    delta = f(x_) - f(x)
    return delta / h_step


# calculate derivatite at x for values of h
f_prime = []
for h_step in h:
    f_prime.append(forward_difference(x_0, h_step))


# print(f"f_prime at {x_0}: {f_prime}")


# Question 1b
def derivative(x):
    x_sq = x ** 2
    return -2 * x * np.exp(-x_sq)


# computing analytic derivative
analytic = derivative(x=0.5)

# print results as required
print("#### Question 1b ####")
for i in range(len(f_prime)):
    # error
    error = abs(f_prime[i] - analytic)

    left_aligned = f"[h={h[i]}]"
    center = f"value:{f_prime[i]}"
    right_aligned = f"error:{error:.16f}"
    print(f"{left_aligned:<12}{center:^30}{right_aligned:>22}")

# Question 1c
# calculating error for plotting
error_arr = []
for i in range(len(f_prime)):
    error_arr.append(abs(f_prime[i] - analytic))

# plotting figure
plt.figure()
indexing = len(h) // 2
plt.scatter(h[indexing], error_arr[indexing], color='red', zorder=2)
plt.plot(h, error_arr, zorder=1, label = "Forward Difference")
plt.annotate(f"{(h[indexing], round(error_arr[indexing],13))}",
             (h[indexing], error_arr[indexing]))
plt.xscale('log')
plt.yscale('log')
plt.title("Logarithmic graph representing error as a function of step size")
plt.xlabel("Step size")
plt.ylabel("Error from differentiating")

# Answering questions
print("\nQuestion 1c")
print("The shape of the curve can be given by a quadratic graph."
      "\nRounding error occurs at very small values of h,"
      "\nwhile trauncation errors occur at very large values of h.")


# Question 1d
def central_difference(x, h_step):
    x_ = x + h_step / 2
    _x = x - h_step / 2
    delta = f(x_) - f(_x)
    return delta / h_step

# calculate derivatite at x for values of h using central difference
f_prime = []
for h_step in h:
    f_prime.append(central_difference(x_0, h_step))

# calculating error for plotting
error_arr = []
for i in range(len(f_prime)):
    error_arr.append(abs(f_prime[i] - analytic))

# plotting figure
indexing = len(h) // 2 + 2
plt.scatter(h[indexing], error_arr[indexing], color='green', zorder=2)
plt.plot(h, error_arr, zorder=1, label = "Central Difference")
plt.annotate(f"{(h[indexing], round(error_arr[indexing],13))}",
             (h[indexing], error_arr[indexing]))
plt.legend()
plt.grid()
plt.savefig("fig1.pdf")
