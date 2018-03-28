import matplotlib.pyplot as plt    
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.colors import ListedColormap, BoundaryNorm    
import numpy as np
from scipy.integrate import ode
n = 500
dt = 10./(n)
rho = 28.
sigma = 10.
beta = 8./3.
x = np.zeros(n)
y = np.zeros(n)
z = np.zeros(n)

u0 = [10,10,10]

def derivative(x,y,z,dt):
    vx = np.zeros(len(x)-1)
    vy = np.zeros(len(x)-1)
    vz = np.zeros(len(x)-1)

    for i in range(len(x)-1):
        vx[i]= (x[i+1] - x[i])/dt
        vy[i]= (y[i+1] - y[i])/dt
        vz[i]= (z[i+1] - z[i])/dt
    return vx,vy,vz


def f(t,u,sigma,rho,beta):
    x,y,z = u
    return [sigma*(y-x), x*(rho-z) - y, x*y-beta*z]



r1 = ode(f).set_integrator('dopri5')
r1.set_initial_value(u0, 0).set_f_params(sigma,rho,beta)
i=0
while  i < n: #r1.successful()  
    x[i],y[i],z[i] = r1.integrate(r1.t+dt)
    #print("%g %g" % (r1.t, r1.y))
    i+=1


#post integration (x, y and z is defined):


vx = sigma*(y-x)
vy = x*(rho-z)-y
vz = x*y-beta*z
x_2 = sigma*(vy - vx) 
y_2 = vx*(rho - z) - x*vz - vy
z_2 = vx*y + x*vy - beta*vz
x_3 = sigma*(y_2 - x_2)
y_3 = x_2*(rho - z) - 2*vx*vz - x*z_2 - y_2
z_3 = x_2*y + 2*vx*vy + x*y_2 - beta*z_2

curvature = np.zeros(n) #empty arrays
torsion   = np.zeros(n)
for i in range(n):
    #first we define the elements in the equations
    d  = np.asarray([ vx[i] ,vy[i] ,vz[i]]) #r dot
    dd = np.asarray([x_2[i],y_2[i],z_2[i]]) #r dotdot
    ddd = np.asarray([x_3[i],y_3[i],z_3[i]])#r dotdotdot
    dxdd = np.cross(d,dd)                   #rdot x rdotdot
    dnorm = np.linalg.norm(d)               #||rdot||
    dxddnorm = np.linalg.norm(dxdd)         #||rdot x rdotdot|| 
    
    curvature[i] = dxddnorm/(dnorm*dnorm*dnorm)
    torsion[i] = np.dot(d, np.cross(dd,ddd)) / (dxddnorm*dxddnorm)



def update_line(num, r,line,ax,segments):
    
    line.set_data(r[0:2, :num])
    line.set_3d_properties(r[2, :num])
    #all the lc and ax lines should not be inside here for normal animation
    lc = Line3DCollection(segments[0:(num+1)], cmap=plt.get_cmap('coolwarm'),
                    norm=norm)
    lc.set_array(torsion[0:(num+1)]) 
    lc.set_linewidth(2)

    ax.add_collection3d(lc, zs=r[2,0:num+1], zdir='z')

    return line

r = np.asarray([x,y,z])
fig = plt.figure()
ax = fig.gca(projection='3d')

norm = plt.Normalize(torsion.min(), torsion.max())
points = np.array([x, y, z]).T.reshape(-1, 1, 3)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Create the 3D-line collection object

lines = ax.plot(x[0:1], y[0:1], z[0:1])[0]

#ax.add_collection3d(lc, zs=z, zdir='z')


line_ani = animation.FuncAnimation(fig, update_line, n-1, fargs=(r, lines,ax,segments),
                                   interval=20, blit=False)

ax.set_xlim3d([min(x), max(x)])
ax.set_ylim3d([min(y), max(y)])
ax.set_zlim3d([min(z), max(z)])
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.show()

