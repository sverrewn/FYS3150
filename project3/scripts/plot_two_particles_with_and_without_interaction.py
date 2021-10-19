import matplotlib.pyplot as plt
import numpy as np

infile1 = 'data/two_particles_interaction.txt'
infile2 = 'data/two_particles_no_interaction.txt'

with open(infile1, 'r') as file:
    data = file.readlines()
    p1 = data[1::7]
    p2 = data[4::7]

    p1 = [x.strip() for x in p1]
    p2 = [x.strip() for x in p2]

    x1 = [float(x.split(',')[0]) for x in p1]
    y1 = [float(y.split(',')[1]) for y in p1]

    x2 = [float(x.split(',')[0]) for x in p2]
    y2 = [float(y.split(',')[1]) for y in p2]

with open(infile2, 'r') as file:
    data = file.readlines()
    p1 = data[1::7]
    p2 = data[4::7]

    p1 = [x.strip() for x in p1]
    p2 = [x.strip() for x in p2]

    x3 = [float(x.split(',')[0]) for x in p1]
    y3 = [float(y.split(',')[1]) for y in p1]

    x4 = [float(x.split(',')[0]) for x in p2]
    y4 = [float(y.split(',')[1]) for y in p2]    

plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)

plt.figure()

plt.plot(x1, y1, label='p1 interact')
plt.plot(x2, y2, label='p2 interact')

plt.plot(x3, y3, label='p1 no interact')
plt.plot(x4, y4, label='p2 no interact')
plt.xlabel('x')
plt.ylabel('y')

plt.grid()
plt.legend()
plt.savefig('figs/two_particles_with_and_without_interaction.pdf')