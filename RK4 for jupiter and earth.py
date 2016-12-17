# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 16:48:45 2016

@author: xh810
"""

import numpy as np
from math import pi
import matplotlib.pyplot as plt

npoints=10000
dt=0.001
Ms=2.0e30 #mass of sun in kg
Me=6.0e24 #mass of earth in kg
Mj=1.9e27 #mass of jupiter in kg

def f(r,t):
    x1=r[0]
    y1=r[1]
    vx1=r[2]
    vy1=r[3]
    x2=r[4]
    y2=r[5]
    vx2=r[6]
    vy2=r[7]
    fx1=vx1
    fx2=vx2
    fy1=vy1
    fy2=vy2
    fvx1=-4*pi**2*x1/(x1**2+y1**2)**1.5-4*pi**2*(x1-x2)/((x1-x2)**2+(y1-y2)**2)**1.5
    fvx2=-4*pi**2*x2/(x2**2+y2**2)**1.5-4*pi**2*(x2-x1)/((x1-x2)**2+(y1-y2)**2)**1.5
    fvy1=-4*pi**2*y1/(x1**2+y1**2)**1.5-4*pi**2*(y1-y2)/((x1-x2)**2+(y1-y2)**2)**1.5
    fvy2=-4*pi**2*y2/(x2**2+y2**2)**1.5-4*pi**2*(y2-y1)/((x1-x2)**2+(y1-y2)**2)**1.5
    return np.array([fx1,fy1,fvx1,fvy1,fx2,fy2,fvx2,fvy2],float)
    
x1_initial=1 #initial position of earth in AU
y1_initial=0
vx1_initial=0 #initial velocity of earth in AU/yr
vy1_initial=2*pi

x2_initial=5.2 #initial position of Jupiter in AU; assume at opposition initially
y2_initial=0
vx2_initial=0  #iniital velocity of Jupiter
vy2_initial=2.7549 #2*pi*5.2 AU/11.85 years=2.7549 AU/yr

r=np.array([x1_initial,y1_initial,vx1_initial,vy1_initial,x2_initial,y2_initial,vx2_initial,vy2_initial],float)


tpoints=np.arange(0,npoints+dt,dt)
x1points=[]
y1points=[]
vx1points=[]
vy1points=[]
x2points=[]
y2points=[]
vx2points=[]
vy2points=[]

for t in tpoints:
    x1points.append(r[0])
    y1points.append(r[1])
    vx1points.append(r[2])
    vy1points.append(r[3])
    x2points.append(r[4])
    y2points.append(r[5])
    vx2points.append(r[6])
    vy2points.append(r[7])
    k1=dt*f(r,t)
    k2=dt*f(r+0.5*k1,t+0.5*dt)
    k3=dt*f(r+0.5*k2,t+0.5*dt)
    k4=dt*f(r+k3,t+dt)
    r+=(k1+2*k2+2*k3+k4)/6
x1points.append(r[0])
y1points.append(r[1])
vx1points.append(r[2])
vy1points.append(r[3])
x2points.append(r[4])
y2points.append(r[5])
vx2points.append(r[6])
vy2points.append(r[7])


    
fig=plt.figure()
plt.plot(x1points,y1points,'r.',x2points,y2points,'k.')
plt.xlim(-7,7)
plt.ylim(-7,7)
plt.xlabel('x(AU)')
plt.ylabel('y(AU)')
plt.title('three body simulation-Jupiter*&Earth')
#fig.savefig('jupiter and earth')
plt.show()

