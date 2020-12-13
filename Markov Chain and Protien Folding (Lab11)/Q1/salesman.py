from math import sqrt,exp
from numpy import empty
from random import random,randrange,seed
import matplotlib.pyplot as plt
from matplotlib import rc

font = {'family': 'DejaVu Sans', 'size': 14}  # adjust fonts
rc('font', **font)

seed(10)

N = 25
R = 0.02
Tmax = 10.0
Tmin = 1e-3
tau = 1e5

# Function to calculate the magnitude of a vector
def mag(x):
    return sqrt(x[0]**2+x[1]**2)

# Function to calculate the total length of the tour
def distance():
    s = 0.0
    for i in range(N):
        s += mag(r[i+1]-r[i])
    return s

# Choose N city locations and calculate the initial distance
r = empty([N+1,2],float)
for i in range(N):
    r[i,0] = random()
    r[i,1] = random()
r[N] = r[0]
D = distance()

# # Plot the initial config
# plt.figure()
# plt.plot(r[:,0], r[:,1], '-o')

# Main loop
t = 0
T = Tmax

# smallest = 1
seed_val = 3
seed(seed_val)
while T>Tmin:

    # Cooling
    t += 1
    T = Tmax*exp(-t/tau)

    # # Update the visualization every 100 moves
    # if t%10000==0:
    #     plt.figure()
    #     plt.plot(r[:,0], r[:,1], '-o')
    
    # Choose two cities to swap and make sure they are distinct
    i,j = randrange(1,N),randrange(1,N)
    
    while i==j:
        i,j = randrange(1,N),randrange(1,N)

    # Swap them and calculate the change in distance
    oldD = D
    r[i,0],r[j,0] = r[j,0],r[i,0]
    r[i,1],r[j,1] = r[j,1],r[i,1]
    D = distance()
    deltaD = D - oldD

    # If the move is rejected, swap them back again
    if random()>exp(-deltaD/T):
        r[i,0],r[j,0] = r[j,0],r[i,0]
        r[i,1],r[j,1] = r[j,1],r[i,1]
        D = oldD
        
plt.figure(figsize=(8,6))
plt.plot(r[:,0], r[:,1], '-o')
plt.xlabel('x')
plt.ylabel('y')
plt.title("Final path with seed={}, tau = {}".format(seed_val, tau))
plt.tight_layout()
plt.grid()
plt.savefig('1a_seed_{}_tau_{}.pdf'.format(seed_val, tau))

print("For seed = {} and tau={}, the total distance D = {:.2f}.".format(seed_val, tau, D))

