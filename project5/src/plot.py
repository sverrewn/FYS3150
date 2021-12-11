import numpy as np
import matplotlib.pyplot as plt

import pyarma as pa


U_1 = pa.cx_cube()
U_2 = pa.cx_cube()

U_1.load("data/run1.bin")
U_2.load("data/run2.bin")
U_1 = np.array(U_1).reshape(np.shape(U_1)[0], np.shape(U_1)[1] * np.shape(U_1)[2])
U_2 = np.array(U_2).reshape(np.shape(U_2)[0], np.shape(U_2)[1] * np.shape(U_2)[2])

p_1 = np.sum(np.real(U_1*np.conj(U_1)), axis=1)
p_2 = np.sum(np.real(U_2*np.conj(U_2)), axis=1)



plt.figure()
plt.plot(np.log10(p_1 - 1))
plt.plot(np.log10(p_2 - 1))
plt.show()
