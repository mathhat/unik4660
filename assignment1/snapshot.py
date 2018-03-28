import matplotlib.pyplot as plt    
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.colors import ListedColormap, BoundaryNorm    
import numpy as np
from scipy.integrate import ode
n = 10000
dt = 9./(n)
rho = 14. #28
sigma = 10. #10
beta = 8. #8/3
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
T = np.zeros((3,n))
N = np.zeros((3,n))
B = np.zeros((3,n))

#first we define the elements in the equations
d  = np.asarray([ vx[0] ,vy[0] ,vz[0]]) #r dot
dd = np.asarray([x_2[0],y_2[0],z_2[0]]) #r dotdot
ddd = np.asarray([x_3[0],y_3[0],z_3[0]])#r dotdotdot
dxdd = np.cross(d,dd)                   #rdot x rdotdot
dnorm = np.linalg.norm(d)               #||rdot||
dxddnorm = np.linalg.norm(dxdd)         #||rdot x rdotdot|| 

T[:,0] = d[:]/dnorm
curvature[0] = dxddnorm/(dnorm*dnorm*dnorm)
torsion[0] = np.dot(d, np.cross(dd,ddd)) / (dxddnorm*dxddnorm)

for i in range(1,n):
    #first we define the elements in the equations
    d  = np.asarray([ vx[i] ,vy[i] ,vz[i]]) #r dot
    dd = np.asarray([x_2[i],y_2[i],z_2[i]]) #r dotdot
    ddd = np.asarray([x_3[i],y_3[i],z_3[i]])#r dotdotdot
    dxdd = np.cross(d,dd)                   #rdot x rdotdot
    dnorm = np.linalg.norm(d)               #||rdot||
    dxddnorm = np.linalg.norm(dxdd)         #||rdot x rdotdot|| 
    
    T[:,i] = d[:]/dnorm
    dT = (T[:,i]-T[:,i-1])/dt
    N[:,i] = dT/np.linalg.norm(dT)
    B[:,i] = np.cross(T[:,i],N[:,i])
    curvature[i] = dxddnorm/(dnorm*dnorm*dnorm)
    torsion[i] = np.dot(d, np.cross(dd,ddd)) / (dxddnorm*dxddnorm)
B = 3*B
N=3*N
T=3*T

def update_line(num, r,line,ax):
    
    line[0].set_data(r[0:2, :num])
    line[0].set_3d_properties(r[2, :num])
    line[1].set_data([[x[num],x[num]+B[0,num]],[y[num],y[num]+ B[1,num]]])
    line[1].set_3d_properties([z[num], z[num]+B[2,num]])
    line[2].set_data([[x[num],x[num]+N[0,num]],[y[num],y[num]+ N[1,num]]])
    line[2].set_3d_properties([z[num], z[num]+N[2,num]])
    line[3].set_data([[x[num],x[num]+T[0,num]],[y[num],y[num]+ T[1,num]]])
    line[3].set_3d_properties([z[num], z[num]+T[2,num]])
    #line[2].set_data([r[0:2, num], r[0:2, num]+N[0:2, num]])
    #line[2].set_3d_properties([r[2, num],r[2, num]+N[2, num]])
    #line[3].set_data([r[0:2, num], r[0:2, num]+T[0:2, num]])
    #line[3].set_3d_properties([r[2, num],r[2, num]+T[2, num]])
    
    

    return line




r = np.asarray([x,y,z])
fig = plt.figure()

ax = fig.gca(projection='3d')

ax.set_xlim3d([min(x), max(x)])
ax.set_ylim3d([min(y), max(y)])
ax.set_zlim3d([min(z), max(z)])
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.title(r"Torsion Color Plot, $\rho = 14,\ \sigma = 20,\ \beta = 8$",size=19)
#colors
norm = plt.Normalize(min(torsion), max(torsion)-0.1)
points = np.array([x, y, z]).T.reshape(-1, 1, 3)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

lc = Line3DCollection(segments, cmap=plt.get_cmap('coolwarm'),
                    norm=norm)
lc.set_array(torsion) 
lc.set_linewidth(2)
ax.add_collection3d(lc, zs=z, zdir='z')


'''
line1 = ax.plot(x[0:1], y[0:1], z[0:1])[0]
line2 = ax.plot(x[0:1]+B[0,0:1],y[0:1]+ B[1,0:1], z[0:1]+B[2,0:1],'r',label="B")[0]
line3 = ax.plot(x[0:1]+N[0,0:1],y[0:1]+ N[1,0:1], z[0:1]+N[2,0:1],'b',label="N")[0]
line4 = ax.plot(x[0:1]+T[0,0:1],y[0:1]+ T[1,0:1], z[0:1]+T[2,0:1],'g',label="T")[0]
'''

#lines = [line1,line2,line3,line4]
#ax.plot(x, y, z)

#line_ani = animation.FuncAnimation(fig, update_line, n, fargs=(r, lines,ax),
#                                   interval=5, blit=False)
plt.colorbar(lc)

plt.show()
