import numpy as np
import matplotlib.pyplot as plt

def main():
    infile1 = 'data/single_particle_100us_pos.txt'
    #infile2 = 'data/two_particles_interaction_pos.txt'
    #infile_vel = 'data/two_particles_interaction_vel.txt'
    #infile3 = 'data/two_particles_no_interaction_pos.txt'
    #plot_1_particle(infile1)
    #plt.figure()
    #plot_2_particles(infile2)
    #plot_2_particles(infile3)
    #plt.show()
    #r = read_data(infile2)
    #v = read_data(infile_vel)
    #plot_phase_space(r, v)

    r = read_data(infile1)
    t_end = 100; dt = 1e-3; 
    r_a, t = analytical_solution(500, 300, 110, t_end, dt)
    r_err = np.linalg.norm((r[0] - r_a), axis=1)

    plt.plot(t, r_err)
    plt.show()



def read_data(filename):
    infile = np.loadtxt(filename)
    N = int(np.max(infile.transpose()[0]) + 1)
    index = np.where(infile.transpose()[0] == 0)
    x = np.array(infile[index, 1:])
    for i in range(1, N):
        index = np.where(infile.transpose()[0] == i)
        x = np.append(x, infile[index, 1:], axis=0)
    return x

def plot_1_particle(file):
    r = read_data(file)
    t = np.arange(0, 100.001, 0.001)
    x, y, z = np.transpose(r[0])
    plt.figure()
    plt.plot(t, z)
    plt.savefig('figs/one_particle_100us.pdf')


def plot_2_particles(file1):
    r = read_data(file1)
    r1 = r[0]
    r2 = r[1]

    plt.grid()

    plt.plot(np.transpose(r1)[0], np.transpose(r1)[1])
    plt.plot(np.transpose(r2)[0], np.transpose(r2)[1])


def plot_phase_space(r, v):
    fig, axs = plt.subplots(3)
    for i in range(len(r)):
        x, y, z = np.transpose(r[i])
        v_x, v_y, v_z = np.transpose(v[i])
        axs[0].plot(x, v_x)
        axs[1].plot(y, v_y)
        axs[2].plot(z, v_z)
    plt.show()


def plot_3D(r):
    r = np.transpose(r, (0, 2, 1))

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    for i in range(len(r)):
        plot = [ax.plot3D(r[i, 0, :], r[i, 1, :], r[i, 2, :])]
    plt.show()


def analytical_solution(x_0, z_0, v_0, t_end, dt, B_0 = 9.65e1, V_0 = 9.65e8, m = 40.78, q = 1, d = 1e4):
    t = np.arange(0, t_end + dt, dt)

    omega_0 = q * B_0 / m
    omega_z = np.sqrt(2*q * V_0 / (m*d**2))

    omega_plus = (omega_0 + np.sqrt(omega_0**2 - 2*omega_z**2)) / 2
    omega_minus = (omega_0 - np.sqrt(omega_0**2 - 2*omega_z**2)) / 2

    A_plus = (v_0 + omega_minus*x_0) / (omega_minus - omega_plus)
    A_minus = -(v_0 + omega_plus*x_0) / (omega_minus - omega_plus)

    f = A_plus*np.exp(-1j*omega_plus*t) + A_minus*np.exp(-1j*omega_minus*t)
    return np.transpose(np.array([np.real(f), np.imag(f), z_0*np.cos(omega_z*t)])), t

main()