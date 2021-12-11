import numpy as np
import matplotlib.pyplot as plt

import pyarma as pa

def deviation(filename):
    U = pa.cx_cube()
    U.load(filename)
    U = np.array(U).reshape(np.shape(U)[0], np.shape(U)[1] * np.shape(U)[2])

    p = np.sum(np.real(U*np.conj(U)), axis=1)
    
    plt.plot(np.log10(p - 1))
"""    
plt.figure()
deviation("data/run1.bin")
deviation("data/run2.bin")
plt.show()
"""
def colmap(filename, t, dt):
    i=[]
    for t in t:
        i.append(int(t/dt))
    
    U = pa.cx_cube()
    U.load(filename)
    U = np.array(U)

    p = np.real(U*np.conj(U))
    for i in i:
        plt.contourf(p[i])
        plt.colorbar()
        plt.show()

dt = 2.5e-5
t = [0, 0.001, 0.002]
colmap("data/run2.bin", t, dt)
