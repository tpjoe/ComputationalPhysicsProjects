import numpy as np


def solveTheta(A, p, q):
    tau = (A[q, q] - A[p, p]) / (2 * A[p, q])
    if tau <= 0:
        t = 1.0 / (tau + np.sqrt(1.0 + tau * tau))
    else:
        t = -1.0 / (-tau + np.sqrt(1.0 + tau * tau))
    c = 1/np.sqrt(1+t*t)
    s = c*t
    # theta = np.arctan(t)
    theta = [s, c, t]
    return theta


def main():
    A = np.matrix('1 2 0; 2 1 2; 0 2 1')
    p = 1
    q = 2
    theta = solveTheta(A, p, q)
    print 'theta =', theta
    return


if __name__ == '__main__':
    main()
