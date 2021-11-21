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

    print(data[1][0:half])

    plt.figure()
    plt.plot(np.log10(data[1][0:half]), data[5][0:half])
    plt.plot(np.log10(data[1][half:]), data[5][half:])
    plt.show()

    plt.figure()
    plt.plot(np.log10(data[1][0:half]), data[9][0:half])
    plt.plot(np.log10(data[1][half:]), data[9][half:])
    plt.show()


if __name__ == '__main__':
    #plot()
    pass