
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import time


# Gaussian Elimination
fig1 = plt.figure()
fig2 = plt.figure()

# Define the test case
testcase = (10, 100, 1000, 10000, 200000, 1000000)
print 'Gaussian Elimination'

# Loop for different number of points
for m in testcase:

    # define variables
    n = m
    d1 = 2
    e1 = -1
    h = 1 / (float(n+1) - 1)
    d = d1 * np.ones((1, n-1))  # Diagonal
    e = e1 * np.ones((1, n-2))  # Diagonal-1
    x = (np.arange(n+1.0-2.0)) / n + h  # x
    f = (h ** 2.0) * 100.0 * np.exp(-10.0 * x)  # fi
    u = np.zeros((1, n+1))

    start = time.clock()                   # Start counting time
    # Forward Substitution
    for i in xrange(1, n-1):
        d[0, i] = d[0, i] - e[0, i - 1] / d[0, i - 1] * e[0, i - 1]
        f[i] = f[i] - e[0, i - 1] / d[0, i - 1] * f[i - 1]
    # Backward Substitution
    u[0, n-2] = f[n-2] / d[0, n-2]
    for i in range(n-3, -1, -1):
        u[0, i] = (f[i] - e[0, i] * u[0, i + 1]) / d[0, i]
    finish = time.clock()                  # End counting time
    # Exact solution and error
    u = u[0,range(0,n-1)]                    # Delete corresponding u
    u_exact = 1-(1-np.exp(-10.0))*x-np.exp(-10.0*x)
    e = np.log10(abs(u-u_exact)/u_exact)      # error
    e1max = max(e)
    print 'The maximum error from Gaussian elimination using', m, 'step is', e1max

    # Plot
    plt.figure(1)
    plt.plot(x, u, 'o', label=m)
    plt.figure(2)
    plt.plot(x, e, label=m)
    plt.legend()

    print 'time for Gaussian Elimination is ', finish - start, 's'

plt.figure(1)
plt.title('Calculation results')
plt.plot(x, u_exact, label='exact')
plt.legend()
plt.figure(2)
plt.title('Relative error')



plt.figure(3)
plt.plot(x, u_exact, label='exact')
plt.legend()
plt.show()