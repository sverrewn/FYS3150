import matplotlib.pyplot as plt
import numpy as np

def plot_solve_diff(N):
    
    eigenvectors = []
    
    data = np.genfromtxt(f"data/discretization_n{N}.csv", comments="#", delimiter=",")
        
    eigenvals = data[:, 0]
    eigenvectors = data[:, 1:]
    
    
        
    plt.plot(eigenvectors[0])
    plt.plot(eigenvectors[1])
    plt.plot(eigenvectors[2])
    
    plt.legend([eigenvals[0], eigenvals[1], eigenvals[2]])
    
    plt.title(f"Discretization of x_hat with n = {N}")
    plt.xlabel("Positions")
    plt.ylabel("Vector elements")
    
    plt.savefig(f"plots/discretization_n{N}.pdf")
    
if __name__ == "__main__":
    plot_solve_diff(10)
    plot_solve_diff(100)
    
    