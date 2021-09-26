import matplotlib.pyplot as plt
import numpy as np

def plot_solve_diff(N, anal=False):
    
    eigenvectors = []
    
    data = np.genfromtxt(f"data/discretization_n{N}.csv", comments="#", delimiter=",")
    
    anal_data = np.genfromtxt(f"data/anal_n{N}.csv", comments="#", delimiter=",")    
        
    eigenvals = data[:, 0]
    eigenvectors = data[:, 1:]
    
    anal_eigenvals = anal_data[:, 0]
    anal_eigenvectors = anal_data[:, 1:]
    
    np.pad(eigenvals, (1, 1), 'constant') # Pad 1 left, 1 right
    np.pad(eigenvectors, (1, 1), 'constant')
    
    np.pad(anal_eigenvals, (1, 1), 'constant')
    np.pad(anal_eigenvectors, (1, 1), 'constant')
    
    
    
    x = range(N) 
    
    for i in x:
        i = i * 1/(N+1)
        
    plt.plot(x, eigenvectors[0], 'o-')
    plt.plot(x, eigenvectors[1], 'o-')
    plt.plot(x, eigenvectors[2], 'o-')
    
    plt.plot(x, np.sign(eigenvectors[0] @ anal_eigenvectors[0]) * anal_eigenvectors[0])
    plt.plot(x, np.sign(eigenvectors[1] @ anal_eigenvectors[1]) * anal_eigenvectors[1])
    plt.plot(x, np.sign(eigenvectors[2] @ anal_eigenvectors[2]) * anal_eigenvectors[2])
    
    plt.legend([f'Disc: {eigenvals[0]}', f'Disc: {eigenvals[1]}', f'Disc: {eigenvals[2]}', f'Anal: {anal_eigenvals[0]}', f'Anal: {anal_eigenvals[1]}', f'Anal: {anal_eigenvals[2]}'], loc='upper right')
    
    plt.title(f"Discretization of x_hat with n = {N}")
    plt.xlabel("Positions")
    plt.ylabel("Vector elements")
    
    plt.savefig(f"plots/discretization_n{N}.pdf")
    plt.clf()
    
if __name__ == "__main__":
    plot_solve_diff(10)
    plot_solve_diff(100)
    
    