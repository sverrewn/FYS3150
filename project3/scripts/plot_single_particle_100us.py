import matplotlib.pyplot as plt
import numpy as np

infile = 'data/single_particle_100us.txt'

x = np.arange(0, 100.001, 0.001)
with open(infile, 'r') as file:
    temp = file.readlines()
    temp = temp[1::4]
    temp = [x.strip() for x in temp]
    y = [float(x.split(',')[2]) for x in temp]

plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.figure()
plt.plot(x, y, label='z')

plt.xlabel('t')
plt.ylabel('z')

plt.grid()
plt.legend()
plt.savefig('figs/one_particle_100us.pdf')