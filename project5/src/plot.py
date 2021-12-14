import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pyarma as pa

import time

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
        plt.figure()
        plt.plot(np.log10(p - 1))
        plt.show()

    def plot_wave(self, t, dt):
        i = []
        for t in t:
            i.append(int(t/dt))
        
        for i in i:
            plt.figure()
            plt.contourf(self.p[i])
            plt.colorbar()
            plt.contourf(self.V)
            plt.show()
    
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
        i = []
        for t in t:
            i.append(int(t/dt))
        
        for i in i:
            plt.figure()
            plt.contourf(np.real(self.U[i]))
            plt.colorbar()
            plt.contourf(self.V)
            plt.show()
        
    def plot_imag(self, t, dt):
        i = []
        for t in t:
            i.append(int(t/dt))
        
        for i in i:
            plt.figure()
            plt.contourf(np.imag(self.U[i]))
            plt.colorbar()
            plt.contourf(self.V)
            plt.show()
        
    def plot_measurement(self, t, dt, x = 0.8, h = 0.005):
        i = int(t/dt)
        p = self.p[i].T[int(x/h)]
        plt.figure()
        plt.plot(p / np.sum(p))
        plt.show()

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
        

dt = 2.5e-5
t = [0, 0.001, 0.002]


w = Wave("data/run2.bin")
#w.plot_measurement(0.002, dt)
w.animate_measurement(scale = 3)
