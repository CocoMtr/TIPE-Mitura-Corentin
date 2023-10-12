from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from numpy import *

g = 9.81
m = 55*e-3
k = 0.72
alpha = 1.8*1e-1
lambd = 5

x0, y0, z0 = 0, 0, 0
vx0, vy0, vz0 = 27.36, 19.35, 4.8
wx, wy, wz = 0, 0, 1000

def Euler(a, b, n):
    t = linspace(a, b, n)
    x = [x0] ; y = [y0] ; z = [z0]
    vx = [vx0] ; vy = [vy0] ; vz = [vz0]
    for i in range (n-1):
        dt = t[i+1] - t[i]
        v = sqrt(vx[i]**2 + vy[i]**2 + vz[i]**2)
        vx += [vx[i] - lambd*vx[i]*dt + (alpha/m)*
               (wy*vz[i] - wz*vy[i])*dt]
        vy += [vy[i] - lambd*vy[i]*dt + (alpha/m)*
               (wz*vx[i] - wx*vz[i])*dt]
        vz += [vz[i] - g*dt - lambd*vz[i]*dt +
               (alpha/m)*(wx*vy[i] - wy*vx[i])*dt]
        x += [x[i] + vx[i]*dt]
        y += [y[i] + vy[i]*dt]
        z += [z[i] + vz[i]*dt]
    return(x, y, z)

def Euler_r(a, b, n):
    t = linspace(a, b, n)
    x = [x0] ; y = [y0] ; z = [z0]
    vx = [vx0] ; vy = [vy0] ; vz = [vz0]
    for i in range (n-1):
        dt = t[i+1] - t[i]
        v = sqrt(vx[i]**2 + vy[i]**2 + vz[i]**2)
        vx += [vx[i] - (k/m)*v*vx[i]*dt + (alpha/m)*
               (wy*vz[i] - wz*vy[i])*dt]
        vy += [vy[i] - (k/m)*v*vy[i]*dt + (alpha/m)*
               (wz*vx[i] - wx*vz[i])*dt]
        vz += [vz[i] - g*dt - (k/m)*v*vz[i]*dt +
               (alpha/m)*(wx*vy[i] - wy*vx[i])*dt]
        x += [x[i] + vx[i]*dt]
        y += [y[i] + vy[i]*dt]
        z += [z[i] + vz[i]*dt]
    return(x, y, z)

def tracer_traj(a, n):
    b = a
    while Euler(a, b, n)[2][-1] >= 0:
        b += 1/n
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    ax.plot(Euler(a, b, n)[0], Euler(a, b, n)[1],
            Euler(a, b, n)[2])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

def tracer_traj_r(a, n):
    b = a
    while Euler_r(a, b, n)[2][-1] >= 0:
        b += 1/n
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    ax.plot(Euler_r(a, b, n)[0], Euler_r(a, b, n)[1],
            Euler_r(a, b, n)[2])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    
tracer_traj_r(0, 100)