import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pyarma as pa

class Wave:
    def __init__(self, filename, slits=2):
            self.filename = filename
            self.slits = slits

    @property
    def U(self):
        U = pa.cx_cube()
        U.load(self.filename)
        return np.array(U)

    @property
    def p(self):
        return abs(self.U)**2

    @property
    def V(self):
        V = pa.mat()
        V.load(f"data/potential{self.slits}.bin")
        V = np.array(V)
        return np.where(V != 0, V, np.nan)

    def plot_probalility_deviation(self):
        p = np.sum(self.p.reshape(np.shape(self.p)[0], np.shape(self.p)[1] * np.shape(self.p)[2]), axis=1)
        plt.plot(np.log10(p - 1))

    def save_probalility_deviation(self, outfile):
        plt.figure()
        Wave.plot_probalility_deviation(self)
        plt.savefig(outfile)

    def plot_wave(self, t, dt):
        i = int(t/dt)
        plt.contourf(self.p[i])
        plt.colorbar()
        plt.contourf(self.V)

    def save_wave(self, t, dt, outfile):
        Wave.__save_time_evaluation(self, t, dt, outfile, Wave.plot_wave)

    def animate_wave(self):
        fps = 60
        scale = 1

        p = self.p[0:-1:scale]
        frn = int(len(p)) # number of frames in the animation
        
        fig = plt.figure()

        plot = [plt.contourf(p[0]),]
        cb = plt.colorbar()

        ani = animation.FuncAnimation(fig, Wave.next_frame_wave, frn, fargs=(p, self.V, plot), interval=1000/fps)

        ani.save("animations/slit.mp4", fps=fps)
        plt.close()

    def next_frame_wave(frame_number, p, V, plot):
        plot[0] = [plt.contourf(p[frame_number]),]
        plot[0] = [plt.contourf(V),]

    def plot_real(self, t, dt):
        i = int(t/dt)
        plt.contourf(np.real(self.U[i]))
        plt.colorbar()
        plt.contourf(self.V)
    
    def save_real(self, t, dt, outfile):
        Wave.__save_time_evaluation(self, t, dt, outfile, Wave.plot_real)
        
    def plot_imag(self, t, dt):
        i = int(t/dt)
        plt.contourf(np.imag(self.U[i]))
        plt.colorbar()
        plt.contourf(self.V)
    
    def save_imag(self, t, dt, outfile):
        Wave.__save_time_evaluation(self, t, dt, outfile, Wave.plot_imag)
        
    def plot_measurement(self, t = 0.002, dt = 2.5e-5, x = 0.8, h = 0.005):
        i = int(t/dt)
        p = self.p[i].T[int(x/h)]
        plt.plot(p / np.sum(p))
    
    def save_measurement(self, outfile, t = 0.002, dt = 2.5e-5, x = 0.8, h = 0.005):
        plt.figure()
        Wave.plot_measurement(self, t = t, dt = dt, x = x, h = h)
        plt.savefig(outfile)

    def animate_measurement(self, x = 0.8, h = 0.005, scale = 1):
        p = self.p[0:-1:scale]
        fps = 60
        
        fig = plt.figure()
        ax = plt.axes(xlim=(0, len(p[0, 0])))
        line, = ax.plot([], [], lw=2)
        line.set_data([], [])

        y = np.linspace(0, 200, len(self.p[0,0]))

        ani = animation.FuncAnimation(fig, Wave.next_frame_measurement, fargs=(x, h, p, y, line), frames=len(p), interval=1000/fps)

        ani.save("animations/measurement.mp4", fps=fps)
        plt.close()

    def next_frame_measurement(frame_number, x, h, p, y, line):
        p_i = p[frame_number].T[int(x/h)]
        p_i = p_i / np.sum(p_i)
        line.set_data(y, p_i)
        return line,

    def __save_time_evaluation(self, t, dt, outfile, F):
        if isinstance(t, list) == False:
            t = [t]
            outfile = [outfile]
        print(outfile)
        if len(t) != len(outfile):
            raise Exception("Number of outfiles must coencide with number of plots!")

        for i, t in enumerate(t):
            plt.figure()
            F(self, t, dt)
            plt.savefig(outfile[i])
        

dt = 2.5e-5
t = [0, 0.001, 0.002]


w = Wave("data/run2.bin")
w.save_wave(t, dt, ["1.pdf", "2.pdf", "3.pdf"])


"""Manual
First create an instance with the name of a file and the number of slits used in the computation (defaults to 2).


Overview of relevant methods, call on these to get relevant plots. Remember to use a suitable number of outfile names, theese whould end with ".pdf":
save_probalility_deviation(outfile):
    Saves a plot of the deviation from 1 for the total probability.

save_wave(t, dt, outfile):
    Saves plots of the probability for a given time t, t can be a float or a list of floats.

save_real(t, dt, outfile):
    Saves plots of the real part of U for a given time t, t can be a float or a list of floats.

save_imag(t, dt, outfile):
    Saves plots of the imaginary part of U for a given time t, t can be a float or a list of floats.

save_measurement(outfile, t = 0.002, dt = 2.5e-5, x = 0.8, h = 0.005):
    Saves a plot of the probabilty for finding the particle on different places on the detector screen. Arguments defaults to the given values.



If you find it relevant to plot data from several files in the same plot, simply create instances for the different files and call on the "plot" methods instead of the "save" methods (listed above) then save the figures manually with plt.savefig(outfile).
"""