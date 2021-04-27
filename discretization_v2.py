"""
@author: Adi Pall - 4/24/2021
Discretization V2 - using sparse linalg -> increased speed by 30-fold for 101 pts!
Last Updated on Mon Apr 26 22:44:48 2021
"""
from scipy import sparse
import numpy as np
import time

import matplotlib.pyplot as plt

def f(x, y):
    return 0

def bc(x, y, N):
    lBC = 25*np.ones((1,N)).flatten()
    leftBC = lBC[1:N-1]
    
    rBC = 25*np.ones((1,N)).flatten()
    rightBC = rBC[1:N-1]
    
    tBC = 25*X[0,:]
    topBC = tBC[1:N-1]
    
    bBC = 25*np.ones((1,N)).flatten()
    bottomBC = bBC[1:N-1]
    
    lex = np.zeros(((N-2)**2, 1)).flatten()
    
    # Top & Bot BC more challenging due to the way dx corresponds to i, dy corresponds to j
    for i in range(N-2):
        lex[(N-2)*i] = topBC[i] 
    
    for j in range(N-2):
        lex[(N-2)*(j+1)-1] = bottomBC[j]
    
    k1 = np.zeros((len(leftBC),1))
    k1[0] = 1 # place at start of vector
    leftBC_ext = sparse.kron(k1,leftBC).toarray().flatten()
    
    k2 = np.zeros((len(rightBC),1))
    k2[-1] = 1 # place at end of vector
    rightBC_ext = sparse.kron(k2,rightBC).toarray().flatten()
    
    include_in_b = lex + leftBC_ext + rightBC_ext
    
    return [include_in_b, lBC, tBC, rBC, bBC]
    

def form_A(N, hx, hy):
    
    alpha = hx**2/hy**2 # uniform at first

    main_diag = 2 * (1 + alpha) * np.ones((N - 2, 1)).ravel()
    off_diag = -1 * np.ones((N - 2, 1)).ravel()
    
    a = main_diag.shape[0]

    diagonals = [main_diag, off_diag, off_diag]

    B = sparse.diags(diagonals, [0, -1, 1], shape=(a, a)) # creates D block
    
    C = sparse.diags([-1*np.ones((N+1, 1)).ravel()], [0], shape=(a,a)) # creates -I block
        
    I1 = sparse.eye(N-2) # creates I block
    
    A1 = sparse.kron(I1,B) # creates N-2 x N-2 block of D matrices
    
    I2 = sparse.diags([1*np.ones((N, 1)).ravel(),1*np.ones((N, 1)).ravel()], [-1,1], shape=(N-2,N-2))
    
    A2 = sparse.kron(I2,C) # creates N-2 x N-2 block of -I matrices
    
    mat = A1 + A2 # combine both blocks

    return mat

if __name__ == '__main__':

    N = 10 # points, uniform at first
    
    start = time.perf_counter() # time method
    
    (x0, xf) = (0.0, 3)
    (y0, yf) = (0.0, 3)
    
    hx = (xf - x0)/(N-1)
    hy = (yf - y0)/(N-1)
    
    x1 = np.linspace(x0, xf, N)
    y1 = np.linspace(y0, yf, N)
    
    X, Y = np.meshgrid(x1, y1)
    
    ff = f(X,Y) # from function contribution
    
    fbc = bc(X, Y, N) # from BC contribution
    
    b = ff*(hx**2) + fbc[0] # need to include BC when doing it this way
    
    A = form_A(N, hx, hy) 
    
    V = sparse.linalg.spsolve(A,b) # GS goes here --------------------
    elapsed = time.perf_counter() - start # time method
    print(f"time elapsed: {elapsed}")
    
    V = V.reshape((N-2, N-2)).T
    
    U = np.zeros((N,N))
    
    U[1:N-1, 1:N-1] = V # fills in inner vals
    U[:,0] = fbc[1] # left BC
    U[0,:] = fbc[2] # top BC
    U[:,N-1] = fbc[3] # right BC -> OVERWRITES TOP BC (but this makes no difference)
    U[N-1,:] = fbc[4] # bottom BC
    
    # Visualize
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    # surf = ax.plot_surface(X[1:N-1, 1:N-1], Y[1:N-1, 1:N-1], U[1:N-1, 1:N-1])
    # goal of the above was to not plot boundaries, to try and match v1 plot
    surf = ax.plot_surface(X, Y, U)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('T')
    plt.tight_layout()
    plt.show()
    
