import numpy as np


def solveTheta(A, p, q):
    theta = 0.0
    tau = (A[q, q] - A[p, p]) / (2 * A[p, q])
    if (tau <= 0):
        t = 1.0 / (tau + np.sqrt(1.0 + tau * tau))
    else:
        t = -1.0 / (-tau + np.sqrt(1.0 + tau * tau))
    theta = np.arctan(t)
    return theta


def main():
    
    return 0


if __name__ == '__main__':
    main()
