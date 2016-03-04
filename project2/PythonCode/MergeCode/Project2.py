from GenerateA import *
from SolveTheta import *
from CalculateB import *
import numpy as np
from CalculateEigV import *
import matplotlib.pyplot as plt
import time

n = 1000
rho_max = 5
err = 1E-14
A = construct_A(n, rho_max)
# print(A)
maxValue = 999
EigV = np.ones([n, n])
maxReturn = A_max(A)
maxValue = maxReturn[2]
maxPosition = maxReturn[0:2]

start = time.clock()                   # Start counting time
while np.abs(maxValue) > err:

    theta = solveTheta(A, maxPosition[0], maxPosition[1])
    B = CalculateB(A, theta, maxPosition[0], maxPosition[1])
    EigV = CalculateEigV(EigV, theta, maxPosition[0], maxPosition[1])
    A = B
    maxReturn = A_max(A)
    maxValue = maxReturn[2]
    maxPosition = maxReturn[0:2]

finish = time.clock()
EigenValues = sorted(np.diag(B))

print "Three lowest Eigen values = ", EigenValues[0:3]
print "The time used for calculation is ", finish-start




# print EigV
