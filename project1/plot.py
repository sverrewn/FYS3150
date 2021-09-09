import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print('Wrong amount of arguments. This program expects one (1) filename')
    sys.exit()

with open(sys.argv[1]) as file:
    x_vals = file.readline().strip().split(',')
    y_vals = file.readline().strip().split(',')

x_vals = [float(n) for n in x_vals]
y_vals = [float(n) for n in y_vals]


plt.rc("xtick", labelsize=13)
plt.rc("ytick", labelsize=13)
plt.figure()
plt.title('Test')
plt.plot(x_vals, y_vals)

plt.xlabel('x')
plt.ylabel('u(x)')

plt.grid()
plt.show()