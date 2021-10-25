import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def main():
    #infile1 = 'data/single_particle_100us_pos.txt'
    infile2 = 'data/two_particles_interaction_pos.txt'
    #infile3 = 'data/two_particles_no_interaction_pos.txt'
    #plot_1_particle(infile1)
    #plt.figure()
    #plot_2_particles(infile2)
    #plot_2_particles(infile3)
    #plt.show()
    r = read_data(infile2)
    plot_3D(r)



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


def plot_3D(r):
    r = np.transpose(r, (0, 2, 1))

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    for i in range(len(r)):
        plot = [ax.plot3D(r[i, 0, :], r[i, 1, :], r[i, 2, :])]
    plt.show()



main()