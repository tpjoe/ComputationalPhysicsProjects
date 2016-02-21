import numpy as np


def construct_A(n, rho_max):

    # define A, h (step), rho[i], and rho_min
    rho_min = 0
    rho = np.ones(n)
    h = (rho_max - rho_min)/np.double(n)
    V = np.ones(n)
    d = np.ones(n)
    e = -np.ones(n-1)*1/np.square(h)


    # constructing_A
    for i in range(n):         # to create rho_i (from i = 0 to n-1
        rho[i] = rho_min + i*h
        V[i] = np.square(rho[i])
        d[i] = 2/np.square(h) + V[i]

    A = np.identity(n)*d + np.diag(e,-1) + np.diag(e,1)
    A_max = A[n-1,n-1]

    return A

def A_max(n, rho_max):

    A = construct_A(n, rho_max)
    A = np.abs(A)
    Max = np.argmax(A)
    Max_ij = np.unravel_index(Max,(n,n))

    return Max_ij


A_max(5,1)

