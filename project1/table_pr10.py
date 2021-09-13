import matplotlib.pyplot as plt
import numpy as np
import sys


if len(sys.argv) != 2:
    print('Wrong amount of arguments. This program  one datafile')
    sys.exit(1)

with open(sys.argv[1], 'r') as file:
    lines = int(file.readline())
    vals_gen = []
    vals_spec = []

    for i in range(lines):
        vals_gen.append(file.readline().strip().split(','))

    for i in range(lines):
        vals_spec.append(file.readline().strip().split(','))

vals_gen  = [[int(n[0]), float(n[1])] for n in vals_gen]
vals_spec = [[int(n[0]), float(n[1])] for n in vals_spec]

columns = ['n', 'general (s)', 'special (s)']
cells = [[n[0], n[1], m[1]] for n,m in zip(vals_gen, vals_spec)]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.axis('off')
table = ax.table(cellText=cells,
                 colLabels=columns,
                 loc='center')
table.scale(0.75,2)

plt.savefig('plots/cmp_gen_spec.pdf')
