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

