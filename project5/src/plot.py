import numpy as np
import matplotlib.pyplot as plt

import pyarma as pa

def deviation(filename):
    U = pa.cx_cube()
    U.load(filename)
    U = np.array(U).reshape(np.shape(U)[0], np.shape(U)[1] * np.shape(U)[2])

    p = np.sum(np.real(U*np.conj(U)), axis=1)
    
    plt.plot(np.log10(p - 1))
    
plt.figure()
deviation("data/run1.bin")
deviation("data/run2.bin")
plt.show()

