import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) != 2:
    print('Wrong amount of arguments. This program expects one (1) filename')
    sys.exit()

file = sys.argv[1]
x_vals = np.loadtxt(file, dtype=float, delimiter=',', max_rows=1)
y_vals = np.loadtxt(file, dtype=float, delimiter=',', skiprows=1, max_rows=1)

plt.rc("xtick", labelsize=13)
plt.rc("ytick", labelsize=13)
plt.figure()
plt.title('Poisson exact')
plt.plot(x_vals, y_vals)

plt.xlabel('x')
plt.ylabel('u(x)')

plt.grid()
plt.savefig('plots/poisson_exact.pdf')
