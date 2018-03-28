import matplotlib.pyplot as plt 
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

n = 10000
r = 1.
t0 = 1.
z0 = np.arange(10)
theta = np.linspace(0, np.pi*2.,10000)
x = r*np.cos(theta)
y = r*np.sin(theta)
def z(t,theta):
    return np.exp(theta*t)

fig = plt.figure()

ax = fig.gca(projection='3d')
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("Streamlines")
#ax.plot(np.zeros(n),np.zeros(n),np.linspace(z0[0],z[-1],n),'r--',label="Z-Axis")
#for i in range(len(z0)):
for t in np.linspace(-0.3,0.3,3):
    ax.plot(x,y,z(t,theta),label="t = %2.1f"%t)

ax.legend(loc="best")
plt.show()