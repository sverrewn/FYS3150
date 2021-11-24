from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def f(L,a):
    return a/L

if __name__ == '__main__':
    L = [40, 60, 80, 100]
    Tc = [2.275, 2.28, 2.285, 2.27]
    Tc = [t - 2.269 for t in Tc]

    b = curve_fit(f, L, Tc, p0=2)
    print(f'a = {b[0][0]}')
    