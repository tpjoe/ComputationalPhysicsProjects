from GenerateA import *
from SolveTheta import *
from CalculateB import *
import numpy as np
from CalculateEigV import *
import matplotlib.pyplot as plt

n = 4
rho_max = 1
err = 1E-14
A = construct_A(n, rho_max)
print(A)
maxValue = 999
EigV = np.ones(n)
maxReturn = A_max(A)
maxValue = maxReturn[2]
maxPosition = maxReturn[0:2]

while np.abs(maxValue) > err:

    theta = solveTheta(A, maxPosition[0], maxPosition[1])
    B = CalculateB(A, theta, maxPosition[0], maxPosition[1])
    EigV = CalculateEigV(EigV, theta, maxPosition[0], maxPosition[1])
    # print 'EigV = ', EigV
    A = B
    maxReturn = A_max(A)
    maxValue = maxReturn[2]
    maxPosition = maxReturn[0:2]

print B

Bdiag = np.diag(B)
plt.plot(EigV)
plt.show()

Bdiag = np.diag(B)
E = Bdiag*EigV

print EigV
