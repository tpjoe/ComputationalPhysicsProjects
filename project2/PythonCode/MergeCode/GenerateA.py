import numpy as np
from copy import copy

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

    A = np.identity(n)*d + np.diag(e, -1) + np.diag(e, 1)
    A_max = A[n-1, n-1]
    EigV = np.ones(n)

    return A

def A_max(A):

    n = int(np.sqrt(A.size))
    A_copy = copy(A)
    np.fill_diagonal(A_copy,0)
    A_copy = np.abs(A_copy)
    Max = np.argmax(A_copy)
    Max_ij = np.unravel_index(Max, (n, n))
    i = Max_ij[0]
    j = Max_ij[1]
    Max = [i,j, A[i,j]]
#    print A_copy
#    print Max_ij
    return Max

def main():
    print construct_A(5,1)
    A_max(construct_A(5,1))

if __name__ == "__main__":
    main()