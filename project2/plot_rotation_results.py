import matplotlib.pyplot as plt
import numpy as np


data = np.genfromtxt("data/rotation_results.csv", delimiter=",", names=["n", "iterations", "converged"])
    
plt.plot(data['n'], data['iterations'])
plt.title("Number of required transformations scale with the matrix size N")
plt.xlabel("Matrix size N")
plt.ylabel("Number of iterations")
plt.savefig("plots/transformations.pdf")