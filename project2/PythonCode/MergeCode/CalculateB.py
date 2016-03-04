import numpy as np
from copy import copy
from SolveTheta import *

def CalculateB(A, theta, p, q):
    B = copy(A)
    sizeA = np.sqrt(A.size)
    s = theta[0]
    c = theta[1]
    B[p, p] = A[p, p]*c*c-2*A[p, q]*c*s+A[q, q]*s*s
    B[q, q] = A[q, q]*c*c+2*A[p, q]*c*s+A[p, p]*s*s
    B[p, q] = 0
    B[q, p] = 0
    for i in range(0, int(sizeA)):
        if (i != p) and (i != q):
            B[i, p] = A[i, p]*c - A[i, q]*s
            B[p, i] = B[i, p]
            B[i, q] = A[i, q]*c + A[i, p]*s
            B[q, i] = B[i, q]

    return B


def main():

    A = np.matrix('1.0 2 3; 2 1 2; 3 2 1')
    p = 0
    q = 2
    theta = solveTheta(A, p, q)
    print 'theta =', theta
    B = CalculateB(A, theta, p, q)
    print B
    return


if __name__ == '__main__':
    main()
