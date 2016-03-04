import numpy as np
from copy import copy
#from SolveTheta import *

def CalculateEigV(EigV, theta, p, q):

    sizeEigV = np.sqrt(EigV.size)
    s = theta[0]
    c = theta[1]
    for i in range(0,int(sizeEigV)):
        if (i!= p & i!=q):
            EigVp = EigV[i, p]
            EigVq = EigV[i, q]
            EigV[i,p] = s*EigVp-c*EigVq
            EigV[i,q] = s*EigVp+c*EigVq
    return EigV

def main():
    print CalculateEigV(np.ones(5),[0.5,0.5], 2, 3)

if __name__ == '__main__':
    main()