import matplotlib.pyplot as plt
import numpy as np
import sys


def abs_error(exact_val, est_val):
    return abs(exact_val - est_val)


# This doesn't work in our boundary points since they're 0
def rel_error(exact_val, est_val):
    return abs((exact_val - est_val) / (exact_val))


def plot_error(err_type):
    plt.rc("xtick", labelsize=13)
    plt.rc("ytick", labelsize=13)
    plt.figure()
    plt.title(f'{err_type} error')

    for e_file, a_file in zip(e_files, a_files):
        x_vals = np.loadtxt(e_file, dtype=float, delimiter=',', max_rows=1)  # The x values should be the same for both
        y_e_vals = np.loadtxt(e_file, dtype=float, delimiter=',', skiprows=1, max_rows=1)
        y_a_vals = np.loadtxt(a_file, dtype=float, delimiter=',', skiprows=1, max_rows=1)
        
        if err_type == 'absolute':
            plt.semilogy(x_vals, abs_error(y_e_vals, y_a_vals), label=f'n={a_file[16:-4]}')
            plt.ylabel(r'$\log_{10}(\Delta_i$')
        elif err_type == 'relative':
            plt.semilogy(x_vals[1:-1], rel_error(y_e_vals[1:-1], y_a_vals[1:-1]), label=f'n={a_file[16:-4]}')
            plt.ylabel(r'$\log_{10}(\epsilon_i)$')
        else:
            print('Unknown error type. terminating')
            sys.exit()
    plt.xlabel('x')
    plt.legend()

    plt.grid()
    plt.savefig(f'plots/{err_type}_error.pdf')


if len(sys.argv) < 2:
    print('Too few arguments. This program takes n data files in the format exact ... approx ...')
    sys.exit(1)

sep_ind = (len(sys.argv) - 1) / 2  # The index for the last exact file
sep_ind = int(sep_ind) + 1  # If called correctly this will always be an integer

e_files = sorted(sys.argv[1:sep_ind])  # Exact
a_files = sorted(sys.argv[sep_ind:])  # Approx

plot_error('absolute')
plot_error('relative')