import numpy as np
import matplotlib.pyplot as plt
from sys import argv


"""
Data structure:
data[0] = Ordered
data[1] = Cycles
data[2] = Temperature
data[3] = Current E
data[4] = Current M
data[5] = <eps>
data[6] = <E^2>
data[7] = <M>
data[8] = <M^2>
data[9] = <|m|>
data[10] = Cv
data[11] = X
"""

quarter = int((len(argv) -1) / 4)

temp = argv[1]
i = 0

while temp[i] != 'T':
    i += 1
T = [float(arg[i+1:i+5]) for arg in argv[1:quarter + 1]]


data = np.array([np.loadtxt(argv[1])])
i = 2
while True:
    try:
        tmp = np.loadtxt(argv[i])
        data = np.append(data, [tmp], axis=0)
        i += 1

    except:
        data = data.T
        break

def plot():
    plt.figure()
    plt.plot(T, data[5][0:quarter], label=r'L=40')
    plt.plot(T, data[5][quarter:2*quarter], label=r'L=60')
    plt.plot(T, data[5][2*quarter:3*quarter], label=r'L=80')
    plt.plot(T, data[5][3*quarter:4*quarter], label=r'L=100')
    plt.xlabel(r'$T \cdot 1000$')
    plt.ylabel(r'$\langle\epsilon\rangle$')
    plt.legend()
    plt.grid()
    plt.savefig('figs/<e>_differing_L.pdf')

    plt.figure()
    plt.plot(T, data[9][0:quarter], label=r'L=40')
    plt.plot(T, data[9][quarter:2*quarter], label=r'L=60')
    plt.plot(T, data[9][2*quarter:3*quarter], label=r'L=80')
    plt.plot(T, data[9][3*quarter:4*quarter], label=r'L=100')
    plt.xlabel(r'$T \cdot 1000$')
    plt.ylabel(r'$\langle|m|\rangle$')
    plt.legend()
    plt.grid()
    plt.savefig('figs/<|m|>_differing_L.pdf')

    plt.figure()
    plt.plot(T, data[10][0:quarter], label=r'L=40')
    plt.plot(T, data[10][quarter:2*quarter], label=r'L=60')
    plt.plot(T, data[10][2*quarter:3*quarter], label=r'L=80')
    plt.plot(T, data[10][3*quarter:4*quarter], label=r'L=100')
    plt.xlabel(r'$T \cdot 1000$')
    plt.ylabel(r'$C_v$')
    plt.legend()
    plt.grid()
    plt.savefig('figs/Cv_differing_L.pdf')

    plt.figure()
    plt.plot(T, data[11][0:quarter], label=r'L=40')
    plt.plot(T, data[11][quarter:2*quarter], label=r'L=60')
    plt.plot(T, data[11][2*quarter:3*quarter], label=r'L=80')
    plt.plot(T, data[11][3*quarter:4*quarter], label=r'L=100')
    plt.xlabel(r'$T \cdot 1000$')
    plt.ylabel(r'$\chi$')
    plt.legend()
    plt.grid()
    plt.savefig('figs/X_differing_L.pdf')

if __name__ == '__main__':
    plot()