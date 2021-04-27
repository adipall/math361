# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 10:46:14 2021
Discretization V1
@author: Adi Pall
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
import time

def form_A(xr,yt,M=5,N=5):
    """
    Inputs:
        xr: Length (L) of rectangle
        yt: Height (W) of rectange
        M: # of steps in x-dir (default = 5)
        N: # of steps in y-dir (default = 5)
    """
    f = lambda x,y: 0
    g1 = lambda x: 25
    g2 = lambda x: 25*x
    g3 = lambda x: 25
    g4 = lambda x: 25
    # set bottom and left ends of domain at origin to keep it simple
    yb = 0; xl = 0;
    m = M+1; n = N+1;
    x = np.linspace(xl, xr, m)
    y = np.linspace(yb, yt, n)
    h = x[1]-x[0] # step size in x
    k = y[1]-y[0] # step size in y
    h2 = h**2; k2 = k**2;
    sz = m*n
    A = np.zeros([sz,sz])
    b = np.zeros(sz)
    # interior pts
    for i in range(1,M):
        for j in range(1,N):
            A[i+(j)*m,i-1+(j)*m]=1/h2
            A[i+(j)*m,i+1+(j)*m]=1/h2
            A[i+(j)*m,i+(j)*m]=-2/h2-2/k2
            A[i+(j)*m,i+(j-1)*m]=1/k2
            A[i+(j)*m,i+(j+1)*m]=1/k2
            b[i+(j)*m]=f(x[i],y[j])
    # bottom and top boundaries       
    for i in range(0,m):
        j = 0
        A[i+(j)*m,i+(j)*m] = 1
        b[i+(j)*m]= g1(x[i])
        j = N
        A[i+(j)*m,i+(j)*m] = 1
        b[i+(j)*m]= g2(x[i])
    # left and right boundaries
    for j in range(1,N):
        i = 0
        A[i+(j)*m,i+(j)*m] = 1
        b[i+(j)*m]= g3(y[j])        
        i = M
        A[i+(j)*m,i+(j)*m] = 1
        b[i+(j)*m]= g4(y[j])
    return A,b,x,y

if __name__ == '__main__':
    L = 3
    Nx_s = 9 # steps
    Ny_s = 9
    
    start = time.perf_counter()
    
    A,b,x,y = form_A(L,L,Nx_s,Ny_s)
    
    A_sparse = sparse.csr_matrix(A)
    b_sparse = (sparse.csr_matrix(b)).T
    v = sparse.linalg.spsolve(A_sparse, b_sparse)
    elapsed = time.perf_counter() - start
    print(f"time elapsed: {elapsed}")
    w = np.reshape(v,(Nx_s+1,Ny_s+1))
    
    # Visualize
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(x,y,w)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('T')
    plt.tight_layout()
    plt.show()