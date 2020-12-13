import numpy as np
from time import perf_counter
import matplotlib.pyplot as plt
from SolveLinear import GaussElim, PartialPivot
from numpy.random import rand
from numpy.linalg import solve
from tqdm import tqdm


min_size = 5 # min size of the matrix
max_size = 100 # max size of the matrix
# to store the size of matrix
n = []
# to store Gauss elimination timing
gauss_elim_timing = []
# to store Gauss elimination error
gauss_elim_err = []
# to store partial_pivot timing
partial_pivot_timing = []
# to store partial pivot err
partial_pivot_err = []
# to store LUdecomposition timing
LU_decomp_timing = []
# to store LUdecomposition err
LU_decomp_err = []

# run the experiment for matrix of sizes min_size to max_size. tqdm estimates  
# the time for the loop to complete
for N in tqdm(range(min_size, max_size + 1)):
    n.append(N)
    # initialize random A and v
    v = rand(N)
    A = rand(N,N)
    
    # Gauss Elimination
    # start the counter for gauss elimination
    s1 = perf_counter()
    # perform Gauss elimination
    x1 = GaussElim(A, v)
    # Get timing for GaussElim in ms
    gauss_elim_timing.append((perf_counter() - s1) * 1000) 
    # get v by taking dot product of A and x1
    v_sol1 = np.dot(A, x1)
    # append the error which is given by mean(abs(v-v_sol1)) to gauss_elim_err
    gauss_elim_err.append(np.mean(np.abs(v - v_sol1)))
    
    # Partial pivoting
    # start the counter for partial pivoting
    s2 = perf_counter()
    # perform partial pivoting
    x2 = PartialPivot(A, v)
    # Get timing for partial pivoting in ms
    partial_pivot_timing.append((perf_counter() - s2) * 1000) 
    # get v by taking dot product of A and x2
    v_sol2 = np.dot(A, x2)
    # append the error which is given by mean(abs(v-v_sol2)) to partial_pivot_err
    partial_pivot_err.append(np.mean(np.abs(v - v_sol2)))
    
    # LU decomposition
    # start the counter for LU decomposition
    s3 = perf_counter()
    # perform LU decomposition
    x3 = solve(A, v)
    # Get timing for LU decomposition in ms
    LU_decomp_timing.append((perf_counter() - s3) * 1000) 
    # get v by taking dot product of A and x3
    v_sol3 = np.dot(A, x3)
    # append the error which is given by mean(abs(v-v_sol3)) to LU_decomp_err
    LU_decomp_err.append(np.mean(np.abs(v - v_sol3)))


# plot the timing for all three methods
plt.figure(figsize=(8,6))
plt.loglog(n, gauss_elim_timing, label="Gaussian elimination")
plt.loglog(n, partial_pivot_timing, label="Partial pivoting")
plt.loglog(n, LU_decomp_timing, label="LU Decomposition")
plt.legend()
plt.grid()
plt.xlabel("Size of matrix (N)")
plt.ylabel("Time (ms)")    
plt.title("Timing for various methods to solve the random linear system Ax=v")
plt.savefig("Q1bi.pdf")


# plot the error for all three methods
plt.figure(figsize=(8,6))
plt.loglog(n, gauss_elim_err, label="Gaussian elimination")
plt.loglog(n, partial_pivot_err, label="Partial pivoting")
plt.loglog(n, LU_decomp_err, label="LU Decomposition")
plt.legend()
plt.grid()
plt.xlabel("Size of matrix (N)")
plt.ylabel("Error")   
plt.title("Error in solving the random linear system Ax=v using various methods")
plt.savefig("Q1bii.pdf")
    