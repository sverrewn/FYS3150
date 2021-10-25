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

plt.grid()
plt.figure()

plt.plot(e1, r1, label='f = 0.1')
plt.plot(e2, r2, label='f = 0.4')
plt.plot(e3, r3, label='f = 0.7')

plt.xlabel(r'$\omega_V$', fontsize=14)
plt.ylabel('Remaining particles', fontsize=14)
plt.legend()

plt.savefig('figs/remaining_particles.pdf')

with open(infile2, 'r') as file:
    data = file.readlines()

row1 = data[0:75]
row1 = [x.split() for x in row1]
for i in range(len(row1)):
    row1[i] = [float(row1[i][0]), float(row1[i][1]), int(row1[i][2])]
row1.sort(key=lambda x: x[1])


row2 = data[75:150]
row2 = [x.split() for x in row2]
for i in range(len(row2)):
    row2[i] = [float(row2[i][0]), float(row2[i][1]), int(row2[i][2])]
row2.sort(key=lambda x: x[1])

row3 = data[150:225]
row3 = [x.split() for x in row3]
for i in range(len(row3)):
    row3[i] = [float(row3[i][0]), float(row3[i][1]), int(row3[i][2])]
row3.sort(key=lambda x: x[1])

e1 = []; e2 = []; e3 = []; r1 = []; r2 = []; r3 = [];

for row in row1:
    e1.append(row[1])
    r1.append(row[2])
for row in row2:
    e2.append(row[1])
    r2.append(row[2])
for row in row3:
    e3.append(row[1])
    r3.append(row[2])

plt.figure()

plt.plot(e1, r1, label='f = 0.1')
plt.plot(e2, r2, label='f = 0.4')
plt.plot(e3, r3, label='f = 0.7')

plt.xlabel(r'$\omega_V$', fontsize=14)
plt.ylabel('Remaining particles')
plt.legend()

plt.grid()
plt.savefig('figs/remaining_particles_narrow_interaction.pdf')

with open(infile3, 'r') as file:
    data = file.readlines()

row1 = data[0:75]
row1 = [x.split() for x in row1]
for i in range(len(row1)):
    row1[i] = [float(row1[i][0]), float(row1[i][1]), int(row1[i][2])]
row1.sort(key=lambda x: x[1])


row2 = data[75:150]
row2 = [x.split() for x in row2]
for i in range(len(row2)):
    row2[i] = [float(row2[i][0]), float(row2[i][1]), int(row2[i][2])]
row2.sort(key=lambda x: x[1])

row3 = data[150:225]
row3 = [x.split() for x in row3]
for i in range(len(row3)):
    row3[i] = [float(row3[i][0]), float(row3[i][1]), int(row3[i][2])]
row3.sort(key=lambda x: x[1])

e1 = []; e2 = []; e3 = []; r1 = []; r2 = []; r3 = [];

for row in row1:
    e1.append(row[1])
    r1.append(row[2])
for row in row2:
    e2.append(row[1])
    r2.append(row[2])
for row in row3:
    e3.append(row[1])
    r3.append(row[2])

plt.figure()

plt.plot(e1, r1, label='f = 0.1')
plt.plot(e2, r2, label='f = 0.4')
plt.plot(e3, r3, label='f = 0.7')

plt.xlabel(r'$\omega_V$', fontsize=14)
plt.ylabel('Remaining particles')
plt.legend()

plt.grid()
plt.savefig('figs/remaining_particles_narrow_no_interaction.pdf')