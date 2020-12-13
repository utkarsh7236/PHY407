import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter
from tqdm import tqdm


def time_matrix_mult(A, B, n):
    """
    A : n by n constant matrix
    B : n by n constant matrix
    Return time taken to multiply n by n matrices A and B.
    """
    C = np.zeros([n, n], float)
    
    s1 = perf_counter()
    for i in range(n):
        for j in range(n):
            for k in range(n):
                C[i, j] += A[i,k] * B[k,j]
    
    t1 = perf_counter() - s1
    
    return t1
            

def experiment():
    """
    Plot time taken to multiply n by n constant matrices using time_matrix_mult
    and np.dot.
    """
    n_start = 2
    n_end = 250
    time_loop = []
    time_numpy = []
    
    for n in tqdm(range(n_start, n_end)):
        A = np.ones([n,n], float) * 3
        B = np.ones([n,n], float) * 3
        time_loop.append(time_matrix_mult(A, B, n) * 1000)
        
        # time matrix multiplication using np.dot
        s = perf_counter()
        C_matrix = np.dot(A, B)
        t_np = (s - perf_counter()) * 1000
        time_numpy.append(t_np)
    
    plt.figure()
    plt.plot(np.arange(n_start, n_end), time_loop, label="loop implementation")
    plt.plot(np.arange(n_start, n_end), time_numpy, label="numpy.dot")
    plt.grid()
    plt.xlabel("N")
    plt.ylabel("Time (ms)")
    plt.title("Time taken for NxN matrix multiplication")
    plt.legend()
    plt.savefig("q3a.pdf")
    
    
    plt.figure()
    plt.plot(np.arange(n_start, n_end)**3, time_loop, label="loop implementation")
    plt.plot(np.arange(n_start, n_end)**3, time_numpy, label="numpy.dot")
    plt.grid()
    plt.xlabel("N^3")
    plt.ylabel("Time (ms)")
    plt.title("Time taken for NxN matrix multiplication")
    plt.legend()
    plt.savefig("q3b.pdf")


if __name__ == "__main__":
    experiment()
