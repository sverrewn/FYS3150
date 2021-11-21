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

data = np.array([np.loadtxt('data/approx_distr/' + argv[1])])
i = 2
while True:
    try:
        tmp = np.loadtxt('data/approx_distr/' + argv[i])
        data = np.append(data, [tmp], axis=0)
        i += 1

    except:
        data = data.T
        break

def plot():
    plt.figure()
    plt.hist(data[5])
    plt.show()

if __name__ == '__main__':
    plot()