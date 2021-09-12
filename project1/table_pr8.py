import matplotlib.pyplot as plt
import numpy as np
import sys

def rel_error(exact_val, est_val):
    return abs((exact_val - est_val) / (exact_val))


if len(sys.argv) < 2:
    print('Too few arguments. This program takes n data files in the format exact ... approx ...')
    sys.exit(1)

sep_ind = (len(sys.argv) - 1) / 2  # The index for the last exact file
sep_ind = int(sep_ind) + 1  # If called correctly this will always be an integer

e_files = sorted(sys.argv[1:sep_ind])  # Exact
a_files = sorted(sys.argv[sep_ind:])  # Approx

rel_max = []

for e_file, a_file in zip(e_files, a_files):
    y_e_vals = np.loadtxt(e_file, dtype=float, delimiter=',', skiprows=1, max_rows=1)
    y_a_vals = np.loadtxt(a_file, dtype=float, delimiter=',', skiprows=1, max_rows=1)

    rel_max.append( [f'{a_file[16:-4]}', np.amax( rel_error(y_e_vals[1:-1], y_a_vals[1:-1]) ) ])

columns = ['n', 'max rel error']

fig = plt.figure()
ax = fig.add_subplot(111)
ax.axis('off')
table = ax.table(cellText=rel_max,
                 colLabels=columns,
                 loc='center')
table.scale(0.75,2)

plt.savefig('plots/biggest_rel_err.pdf')