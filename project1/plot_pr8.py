import matplotlib.pyplot as plt
import numpy as np
import sys

def abs_error(exact_val, est_val):
    abs_err = np.log(abs(exact_val - est_val))

def rel_error(exact_val, est_val):
    rel_err = np.log(abs((exact_val - est_val) / (exact_val)))

def plot_error(err_type, x_vals, y_vals):
    plt.rc("xtick", labelsize=13)
    plt.rc("ytick", labelsize=13)
    plt.figure()
    plt.title(f'{err_type}_error')
    plt.plot(x_vals, y_vals)

    plt.xlabel('x')
    plt.ylabel('u(x)')

    plt.grid()
    plt.savefig(f'plots/{err_type}_error.pdf')

if __name__ == "__main__":
    # execute only if run as a script
    exact_val = np.loadtxt('data/poisson_exact.dat', dtype=float, delimiter=',', max_rows=1)    
    est_val = np.loadtxt('data/general_tridiag_100000.dat', dtype=float, delimiter=',', max_rows=1)
    
    abs_error(exact_val, est_val)
    
    