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

data1 = np.loadtxt('data/approx_distr/' + argv[1]).T
data2 = np.loadtxt('data/approx_distr/' + argv[2]).T


def plot():
    plt.figure()
    plt.hist(data1, 25,
        edgecolor='black',
        linewidth=1
        )
    plt.xlabel(r'$T [J/k_B]$')
    plt.ylabel(r'counts')
    plt.savefig('figs/hist_T1.pdf')

    plt.figure()
    plt.hist(data2, 25,
        edgecolor='black',
        linewidth=1
        )
    plt.xlabel(r'$T [J/k_B]$')
    plt.ylabel(r'counts')
    plt.savefig('figs/hist_T2.pdf')

if __name__ == '__main__':
    plot()