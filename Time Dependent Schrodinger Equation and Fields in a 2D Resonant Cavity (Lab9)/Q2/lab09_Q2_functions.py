"""
author: Aslesha

Code adapted from the handout.
"""

import numpy as np
from dcst_for_q3 import dct, idct, dst, idst


def dXYt2(f, X, Y):
    """ 
    Takes DXT along x, then DYT along y. Use dct if X/Y = 'C' and dst if X/Y='S'. 
    IN: f, the input 2D numpy array
        X, C - cosine or S - sine
        Y, C - cosine or S - sine
    OUT: b, the 2D transformed array 
    """
    M = f.shape[0] # Number of rows
    N = f.shape[1] # Number of columns
    
    a = np.zeros((M, N)) # Intermediate array
    b = np.zeros((M, N)) # Final array
    
    # Take transform along x
    for j in range(N):
        # DXT f[:, j] and set as a[:, j]
        if X == 'C':
            a[:, j] = dct(f[:, j])
        elif X == 'S':
            a[:, j] = dst(f[:, j])
        else:
            raise ValueError("Use either 'S' or 'C' for X" ) 
        
    # Take transform along y
    for i in range(M):
        # DYT a[i, :] and set as b[i, :]
        if Y == 'C':
            b[i, :] = dct(a[i, :])
        elif Y == 'S':
            b[i, :] = dst(a[i, :])
        else:
            raise ValueError("Use either 'S' or 'C' for Y" ) 
    
    return b


def idXYt2(b, X, Y):
    """ 
    Takes iDYT along y, then iDXT along x. Use idct if X/Y = 'C' and idst if X/Y='S'. 
    IN: b, the input 2D numpy array
        X, C - cosine or S - sine
        Y, C - cosine or S - sine
    OUT: f, the 2D inverse-transformed array 
    """
    M = b.shape[0] # Number of rows
    N = b.shape[1] # Number of columns
    
    a = np.zeros((M, N)) # Intermediate array
    f = np.zeros((M, N)) # Final array
    
    # Take inverse transform along y
    for i in range(M):
        # iDYT b[i,:] and set as a[i,:]
        if Y == 'C':
            a[i, :] = idct(b[i,:])
        elif Y == 'S':
            a[i, :] = idst(b[i,:])
        else:
            raise ValueError("Use either 'S' or 'C' for Y" ) 
        
    # Take inverse transform along x
    for j in range(N):
        # iDXT a[:,j] and set as f[:,j]
        if X == 'C':
            f[:, j] = idct(a[:,j])
        elif X == 'S':
            f[:, j] = idst(a[:,j])
        else:
            raise ValueError("Use either 'S' or 'C' for X" ) 
        
    return f



