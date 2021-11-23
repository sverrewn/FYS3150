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

data = np.array([np.loadtxt('data/burn_in/' + argv[1])])
i = 2
while True:
    try:
        tmp = np.loadtxt('data/burn_in/' + argv[i])
        data = np.append(data, [tmp], axis=0)
        i += 1

    except:
        data = data.T
        break

half = int(( len(argv) - 1 ) / 2)

def plot():

    plt.figure()
    plt.plot(np.log10(data[1][0:half]), data[5][0:half], label=r'$T=1.0$')
    plt.plot(np.log10(data[1][half:]), data[5][half:], label=r'$T=2.4$')
    plt.xlabel(r'$10^x$ cycles')
    plt.ylabel(r'$\langle\epsilon\rangle$')
    plt.legend()
    plt.grid()
    plt.savefig('figs/burn_inE.pdf')

    plt.figure()
    plt.plot(np.log10(data[1][0:half]), data[9][0:half], label=r'T=$1.0$')
    plt.plot(np.log10(data[1][half:]), data[9][half:], label=r'T=$2.4$')
    plt.xlabel(r'$10^x$')
    plt.ylabel(r'$\langle|m|\rangle$')
    plt.legend()
    plt.grid()
    plt.savefig('figs/burn_in|m|.pdf')


if __name__ == '__main__':
    plot()
    pass