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

        ani = animation.FuncAnimation(fig, Wave.next_frame, frn, fargs=(p, self.V, plot), interval=1000/fps)

        ani.save("slit.mp4", fps=fps)
        plt.close()

    def next_frame(frame_number, p, V, plot):
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
        

dt = 2.5e-5
t = [0, 0.001, 0.002]


w = Wave("data/run2.bin")
w.animate_wave()