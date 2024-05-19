import matplotlib.pyplot as plt
import numpy as np
import sys

data = np.loadtxt('out.txt')
data_length = len(data)
x = np.array([3*i/data_length for i in range(data_length)])
f = np.sin(5*x)
fprime = 5*np.cos(5*x)

plt.plot(x,data, label ="numerical derivative")
plt.plot(x,f, label = "f(x)")
plt.plot(x,fprime, label ="df/dx")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.savefig(sys.argv[1])
plt.close()