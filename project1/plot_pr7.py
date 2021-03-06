import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) < 2:
    print('Too few arguments. This program takes n data files')
    sys.exit()

e_file = 'poisson_exact_1000.dat'  # the poisson exact with n=1000 should be enough for this
files = sorted([arg for arg in sys.argv[1:]])

plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.figure()
plt.title('Poisson estimate')

for file in files:
    x_vals = np.loadtxt(file, dtype=float, delimiter=',', max_rows=1)
    y_vals = np.loadtxt(file, dtype=float, delimiter=',', skiprows=1, max_rows=1)
    plt.plot(x_vals, y_vals, label=f'n={file[16:-4]}')  # get n from file name

x_vals = np.loadtxt(e_file, dtype=float, delimiter=',', max_rows=1)
y_vals = np.loadtxt(e_file, dtype=float, delimiter=',', skiprows=1, max_rows=1)
plt.plot(x_vals, y_vals, label='exact')

plt.xlabel('x')
plt.ylabel('u(x)')

plt.grid()
plt.legend()
plt.savefig('plots/gen_tri_cmp_exact.pdf')