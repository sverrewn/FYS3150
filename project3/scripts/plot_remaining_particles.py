import matplotlib.pyplot as plt
import numpy as np

infile1 = 'data/particles_remaining_for_f0.1_0.4_0.7.dat'
infile2 = 'data/narrow_freq_interact.dat'
infile3 = 'data/narrow_freq_no_interact.dat'

data = np.loadtxt(infile1)

row1 = data[0:115]
row2 = data[115:230]
row3 = data[230:345]
e1 = []; e2 = []; e3 = []; r1 = []; r2 = []; r3 = [];

for row in row1:
    e1.append(row[0])
    r1.append(row[1])
for row in row2:
    e2.append(row[0])
    r2.append(row[1])
for row in row3:
    e3.append(row[0])
    r3.append(row[1])


plt.figure()

plt.plot(e1, r1, label='f = 0.1')
plt.plot(e2, r2, label='f = 0.4')
plt.plot(e3, r3, label='f = 0.7')

plt.xlabel(r'$\omega_V$')
plt.ylabel('Remaining particles')
plt.legend()

plt.savefig('figs/remaining_particles.pdf')

with open(infile2, 'r') as file:
    data = readlines()

row1 = data[0:115]
row2 = data[115:230]
row3 = data[230:345]
e1 = []; e2 = []; e3 = []; r1 = []; r2 = []; r3 = [];

for row in row1:
    e1.append(row[0])
    r1.append(row[1])
for row in row2:
    e2.append(row[0])
    r2.append(row[1])
for row in row3:
    e3.append(row[0])
    r3.append(row[1])

