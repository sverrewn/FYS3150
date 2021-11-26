from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def f(L,a):
    return a/L

if __name__ == '__main__':
    L = [40, 60, 80, 100]
    Tc = [2.275, 2.28, 2.285, 2.27]
    T = [t - 2.269 for t in Tc]

    b = curve_fit(f, L, T, p0=2)
    a = b[0][0]

    T_approx = [tc - a/l for tc,l in zip(Tc,L)]
    
    for t,l in zip(T_approx, L):
        print(f'L = {l}: {t}')
    