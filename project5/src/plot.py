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

    def save_wave(self, outfile, t = [0.0, 0.001, 0.002], dt = 2.5e-5):
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
    
    def save_real(self, outfile, t = [0.0, 0.001, 0.002], dt = 2.5e-5):
        Wave.__save_time_evaluation(self, t, dt, outfile, Wave.plot_real)
        
    def plot_imag(self, t, dt):
        i = int(t/dt)
        plt.contourf(np.imag(self.U[i]))
        plt.colorbar()
        plt.contourf(self.V)
    
    def save_imag(self, outfile, t = [0.0, 0.001, 0.002], dt = 2.5e-5):
        Wave.__save_time_evaluation(self, t, dt, outfile, Wave.plot_imag)
        
    def plot_measurement(self, t, dt, x, h):
        i = int(t/dt)
        p = self.p[i].T[int(x/h)]
        plt.plot(p / np.sum(p))
    
    def save_measurement(self, outfile, t = 0.002, dt = 2.5e-5, x = 0.8, h = 0.005):
        plt.figure()
        Wave.plot_measurement(self, t, dt, x, h)
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
        

if __name__ == "__main__":
    infile_1 = "data/run1.bin"
    infile_2 = "data/run2.bin"
    infile_3 = "data/run3.bin"
    infile_4 = "data/run4.bin"
    infile_5 = "data/run5.bin"

    w_1 = Wave(infile_1, slits=0)
    w_2 = Wave(infile_2, slits=2)
    w_3 = Wave(infile_3, slits=2)
    w_4 = Wave(infile_4, slits=1)
    w_5 = Wave(infile_5, slits=3)

    w_1.save_probalility_deviation("figures/prob_deviation_0_slits.pdf")
    w_2.save_probalility_deviation("figures/prob_deviation_2_slits.pdf")

    t_list = [0, 0.001, 0.002]
    ofile_particle_prob = []
    ofile_U_real = []
    ofile_U_imag = []
    for t in t_list:
        ofile_particle_prob.append(f"figures/particle_prob_at_time_{t}.pdf")
        ofile_U_real.append(f"figures/U_real_at_time_{t}.pdf")
        ofile_U_imag.append(f"figures/U_imag_at_time_{t}.pdf")
    w_3.save_wave(ofile_particle_prob)
    w_3.save_real(ofile_U_real)
    w_3.save_imag(ofile_U_imag)

    w_3.save_measurement(f"figures/detector_screen_with_2_slits.pdf")
    w_4.save_measurement(f"figures/detector_screen_with_1_slits.pdf")
    w_5.save_measurement(f"figures/detector_screen_with_3_slits.pdf")
