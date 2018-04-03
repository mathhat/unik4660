import matplotlib.pyplot as plt    
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

n = 10000
dt = 20./n
fig = plt.figure()
ax = fig.gca(projection='3d')
rho = 28.
sigma = 10.
beta = 8./3.
x = np.zeros(n)
y = np.zeros(n)
z = np.zeros(n)
x[0] = 10
y[0] = 10
z[0] = 10

for i in range(n-1):
    x[i+1] = sigma*(y[i]-x[i])*dt + x[i]
    y[i+1] = (x[i]*(rho - z[i])-y[i])*dt + y[i]
    z[i+1] = (x[i]*y[i] - beta*z[i])*dt + z[i]

ax.plot(x, y, z, label='Lorentz Attractor')
ax.legend()
plt.xlabel("X")
plt.ylabel("Y")
ax.set_zlabel("Z")
plt.show()