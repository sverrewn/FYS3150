import numpy as np
import matplotlib.pyplot as plt

def main():
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)

    infile1 = 'data/single_particle_100us_pos.txt'
    infile2 = 'data/two_particles_interaction_pos.txt'
    infile3 = 'data/two_particles_no_interaction_pos.txt'
    infile3_vel = 'data/two_particles_no_interaction_pos.txt'
    infile2_vel = 'data/two_particles_interaction_vel.txt'
    
    rel_err_ec_infiles = [
        'data/one_particle__ec_dt1_pos.txt',
        'data/one_particle__ec_dt0.1_pos.txt',
        'data/one_particle__ec_dt0.01_pos.txt',
        'data/one_particle__ec_dt0.001_pos.txt',
        'data/one_particle__ec_dt0.0005_pos.txt'
    ]
    rel_err_rk4_infiles = [
        'data/one_particle_rk4_dt1_pos.txt',
        'data/one_particle_rk4_dt0.1_pos.txt',
        'data/one_particle_rk4_dt0.01_pos.txt',
        'data/one_particle_rk4_dt0.001_pos.txt',
        'data/one_particle_rk4_dt0.0005_pos.txt'
    ]
    plot_1_particle(infile1)
    plot_2_particles(infile2, infile3)
   
    plot_phase_space(infile2, infile2_vel, infile3, infile3_vel)

    plot_3D(infile2, infile3)

    plot_relative_error(rel_err_ec_infiles, rel_err_rk4_infiles)


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

    plt.ylabel("z[$\\mu$m]", fontsize=14)
    plt.xlabel("Time[$\\mu$s]", fontsize=14)
    plt.savefig('figs/one_particle_100us_z.pdf')


def plot_2_particles(file1, file2):
    
    plt.figure()
    r = read_data(file1)
    r1 = r[0]
    r2 = r[1]

    plt.plot(np.transpose(r1)[0], np.transpose(r1)[1], label='p1 interaction')
    plt.plot(np.transpose(r2)[0], np.transpose(r2)[1], label='p2 interaction')

    r = read_data(file2)
    r1 = r[0]
    r2 = r[1]

    plt.plot(np.transpose(r1)[0], np.transpose(r1)[1], label='p1 no interaction')
    plt.plot(np.transpose(r2)[0], np.transpose(r2)[1], label='p2 no interaction')

    plt.ylabel("y[$\\mu$m]", fontsize=14)
    plt.xlabel("x[$\\mu$m]", fontsize=14)
    plt.grid()
    plt.legend()
    plt.savefig('figs/two_particles_100us.pdf')


def plot_phase_space(inf1, inf2, inf3, inf4):
    r = read_data(inf1)
    v = read_data(inf2)

    fig, axs = plt.subplots(3)
    for i in range(len(r)):
        x, y, z = np.transpose(r[i])
        v_x, v_y, v_z = np.transpose(v[i])
        axs[0].plot(x, v_x)
        
        axs[1].plot(y, v_y)
        
        axs[2].plot(z, v_z)
        
    axs[0].set_ylabel("x", fontsize=14)
    axs[0].set_xlabel("$v_x$", fontsize=14)
    axs[1].set_ylabel('y', fontsize=14)
    axs[1].set_xlabel("$v_y$", fontsize=14)
    axs[2].set_ylabel('z', fontsize=14)
    axs[2].set_xlabel("$v_z$", fontsize=14)
    #plt.legend()
    plt.savefig('figs/phase_space_interact.pdf')

    r = read_data(inf3)
    v = read_data(inf4)

    fig, axs = plt.subplots(3)
    for i in range(len(r)):
        x, y, z = np.transpose(r[i])
        v_x, v_y, v_z = np.transpose(v[i])
        axs[0].plot(x, v_x)
        
        axs[1].plot(y, v_y)
        
        axs[2].plot(z, v_z)
        
    axs[0].set_ylabel("x", fontsize=14)
    axs[0].set_xlabel("$v_x$", fontsize=14)
    axs[1].set_ylabel('y', fontsize=14)
    axs[1].set_xlabel("$v_y$", fontsize=14)
    axs[2].set_ylabel('z', fontsize=14)
    axs[2].set_xlabel("$v_z$", fontsize=14)
    #plt.legend()
    plt.savefig('figs/phase_space_no_interact.pdf')



def plot_3D(f1, f2):
    r = read_data(f1)
    r = np.transpose(r, (0, 2, 1))

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    for i in range(len(r)):
        plot = [ax.plot3D(r[i, 0, :], r[i, 1, :], r[i, 2, :])]
    plt.savefig('figs/3d_plot_interaction.pdf')

    r = read_data(f2)
    r = np.transpose(r, (0, 2, 1))

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    for i in range(len(r)):
        plot = [ax.plot3D(r[i, 0, :], r[i, 1, :], r[i, 2, :])]
    plt.savefig('figs/3d_plot_no_interaction.pdf')

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

def plot_relative_error(r_ec, r_rk4):
    delmax_ec = []
    plt.figure()
    t_end = 100
    dt_i = [1, 1e-1, 1e-2, 1e-3, 5e-4]
    for dt, f in zip(dt_i, r_ec):
        r = read_data(f)
        r_a, t = analytical_solution(500, 300, 110, t_end, dt)
        r_err = np.linalg.norm(r[0] - r_a, axis=1)
        r_relerr = r_err / np.linalg.norm(r_a, axis=1)

        plt.plot(t, r_relerr, label="dt = " + str(dt))

        delmax_ec.append(np.max(r_err))
    plt.xlabel("Time[$\\mu$s]", fontsize=14)
    plt.ylabel("relative error", fontsize=14)
    plt.legend()
    plt.savefig('figs/relative_error_ec.pdf')
    
    delmax_rk = []
    plt.figure()
    for dt, f in zip(dt_i, r_rk4):
        r = read_data(f)
        r_a, t = analytical_solution(500, 300, 110, t_end, dt)
        r_err = np.linalg.norm(r[0] - r_a, axis=1)
        r_relerr = r_err/np.linalg.norm(r_a, axis=1)

        plt.plot(t, r_relerr, label="dt = " + str(dt))

        delmax_rk.append(np.max(r_err))
    plt.xlabel("Time[$\\mu$s]", fontsize=14)
    plt.ylabel("relative error", fontsize=14)
    plt.legend()
    plt.savefig('figs/relative_error_rk4.pdf')
    
    convergence_ec = 0
    convergence_rk = 0
    for i in range(1, 5):
        convergence_ec += 1/4 * np.log10(delmax_ec[i]/delmax_ec[i-1]) / np.log(dt_i[i]/dt_i[i-1])
        convergence_rk += 1/4 * np.log10(delmax_rk[i]/delmax_rk[i-1]) / np.log(dt_i[i]/dt_i[i-1])
    print(convergence_ec)
    print(convergence_rk)


main()