from GenerateA import *
from SolveTheta import *
from CalculateB import *
import numpy as np


n = 5
rho_max = 5
err = 1E-5
A = construct_A(n, rho_max)
maxValue = 999

while np.abs(maxValue) > err:

    maxReturn = A_max(A)
    maxValue = maxReturn[2]
    print A
    print maxValue
    maxPosition = maxReturn[0:2]
    theta = solveTheta(A, maxPosition[0], maxPosition[1])

    B = CalculateB(A, theta, maxPosition[0], maxPosition[1])
    A = B
    print 'Position = ', maxPosition
print B
