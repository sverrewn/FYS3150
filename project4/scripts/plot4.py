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

data = np.array([np.loadtxt('data/test/' + argv[1])])
i = 2
while True:
    try:
        tmp = np.loadtxt('data/test/' + argv[i])
        data = np.append(data, [tmp], axis=0)
        i += 1

    except:
        data = data.T
        break


def plot():

    beta = 1
    N = 4

    eps_int = -2*np.sinh(8*beta) / (np.cosh(8*beta) + 3)
    eps_int2 = 4 * np.cosh(8*beta) / (np.cosh(8*beta) + 3)

    m_int = (np.exp(8*beta) + 2) / (2*np.cosh(8*beta) + 6)
    m_int2 = (np.exp(8*beta) + 1) / (2*np.cosh(8*beta) + 6)


    eps = np.ones(len(data[1])) * eps_int
    m = np.ones(len(data[1])) * m_int
    C_v = np.ones(len(data[1])) * (1 / N) * (eps_int2 * N**2 - (eps_int * N)**2) #k_B
    X = np.ones(len(data[1])) * (m_int2 - m_int**2)
    
    # <E>
    plt.figure()
    plt.plot(np.log10(data[1]), data[5])
    plt.plot(np.log10(data[1]), eps)
    plt.savefig('figs/testE.pdf')

    #<|m|>
    plt.figure()
    plt.plot(np.log10(data[1]), data[9])
    plt.plot(np.log10(data[1]), m)
    plt.savefig('figs/testM.pdf')

    #Cv
    plt.figure()
    plt.plot(np.log10(data[1]), data[10])
    plt.plot(np.log10(data[1]), C_v)
    plt.savefig('figs/testCv.pdf')

    #X
    plt.figure()
    plt.plot(np.log10(data[1]), data[11])
    plt.plot(np.log10(data[1]), X)
    plt.savefig('figs/testX.pdf')

def plot2():
    plt.figure()


if __name__ == '__main__':
    plot()
    pass