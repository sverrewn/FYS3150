import numpy as np
import matplotlib.pyplot as plt
from sys import argv


data = [np.loadtxt(argv[1])]
i = 2
while True:
    try:
        tmp = np.loadtxt(argv[i])
        data = np.append(data, [tmp], axis=0)
        i += 1

    except:
        data = data.T
        break

"""
Data structure:
data[0] = Ordered
data[1] = Cycles
data[2] = Temperature
data[3] = Current E
data[4] = Current M
data[5] = <E>
data[6] = <E^2>
data[7] = <M>
data[8] = <M^2>
data[9] = <|M|>
data[10] = Cv
data[11] = X
"""

