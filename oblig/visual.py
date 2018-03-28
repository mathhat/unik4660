import matplotlib.pyplot as plt 
import numpy as np 
from matplotlib import cm
N = 127*4 #width in cpp document
L = N/10  #length of lines defined in cpp document
numlines = (127*4/10)**2
#numlines=int(np.floor(np.sqrt(numlines))**2)
x = np.zeros((N,N))
y = np.zeros((N,N))
st= np.zeros((2,2*L+1))
fig = plt.figure()
'''
with open ("x.txt","r") as xfile:
    for j in range(N):
        for i in range(N):
            x[j,i] = float(xfile.readline())

with open ("y.txt","r") as yfile:
    for j in range(N):
        for i in range(N):
            y[j,i] = float(yfile.readline())
'''
with open ("lines.txt","r") as linefile:
    for i in range(numlines):
        for j in range(2*L+1):
            line = linefile.readline().split()
            st[0,j] = float(line[0])
            st[1,j] = float(line[1])
        plt.plot(st[0],st[1],"b-")

#X, Y = np.meshgrid(np.arange(0, N), np.arange(0, N))
#from mpl_toolkits.mplot3d import Axes3D
#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(X, Y, x*x + y*y -2000, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)
plt.title("swag")
#Q = plt.quiver(X, Y, x, y,cmap=cm.coolwarm)
plt.show()
