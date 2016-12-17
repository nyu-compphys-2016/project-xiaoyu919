# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 23:14:28 2016

@author: xh810
"""

import numpy as np
from math import pi,sqrt
import matplotlib.pyplot as plt

npoints=1200
dt=0.1
Ms=2.0e30 #mass of sun in kg
Mj=1.9e27 #mass of jupiter in kg

#gap 3:1, T=3.95,v=2*pi*2.5/3.95=3.977
x1_initial=2.5 
y1_initial=0
v1x_initial=0 
v1y_initial=3.977 

#gap 5:2, T=4.74,v=2*pi*2.825/4.74=3.745
x2_initial=2.825 
y2_initial=0
v2x_initial=0 
v2y_initial=3.745

#gap 7:3, T=5.08,v=2*pi*2.95/5.08=3.649
x3_initial=2.95 
y3_initial=0
v3x_initial=0 
v3y_initial=3.649

#gap 2:1, T=5.925,v=2*pi*3.276/5.925=3.474
x4_initial=3.276 
y4_initial=0
v4x_initial=0 
v4y_initial=3.474

xj_initial=5.2 #initial position of Jupiter in AU; assume at opposition initially
yj_initial=0
vjx_initial=0  #iniital velocity of Jupiter
vjy_initial=2.7549 #2*pi*5.2 AU/11.85 years=2.7549 AU/yr

#create array to store position and velocity of asteroids
x1=np.zeros(npoints+1,float)
y1=np.zeros(npoints+1,float)
v1x=np.zeros(npoints+1,float)
v1y=np.zeros(npoints+1,float)

x2=np.zeros(npoints+1,float)
y2=np.zeros(npoints+1,float)
v2x=np.zeros(npoints+1,float)
v2y=np.zeros(npoints+1,float)

x3=np.zeros(npoints+1,float)
y3=np.zeros(npoints+1,float)
v3x=np.zeros(npoints+1,float)
v3y=np.zeros(npoints+1,float)

x4=np.zeros(npoints+1,float)
y4=np.zeros(npoints+1,float)
v4x=np.zeros(npoints+1,float)
v4y=np.zeros(npoints+1,float)

#create array to store position and velocity of Jupoter
xj=np.zeros(npoints+1,float)
yj=np.zeros(npoints+1,float)
vjx=np.zeros(npoints+1,float)
vjy=np.zeros(npoints+1,float)

r1=np.zeros(npoints+1,float)
r2=np.zeros(npoints+1,float)
r3=np.zeros(npoints+1,float)
r4=np.zeros(npoints+1,float)
#
rj=np.zeros(npoints+1,float)
#
r1j=np.zeros(npoints+1,float)
r2j=np.zeros(npoints+1,float)
r3j=np.zeros(npoints+1,float)
r4j=np.zeros(npoints+1,float)

#initial position and velocity of asteroids
x1[0]=x1_initial
y1[0]=y1_initial
v1x[0]=v1x_initial
v1y[0]=v1y_initial

x2[0]=x2_initial
y2[0]=y2_initial
v2x[0]=v2x_initial
v2y[0]=v2y_initial

x3[0]=x3_initial
y3[0]=y3_initial
v3x[0]=v3x_initial
v3y[0]=v3y_initial

x4[0]=x4_initial
y4[0]=y4_initial
v4x[0]=v4x_initial
v4y[0]=v4y_initial

#initial position and velocity of Jupiter
xj[0]=xj_initial
yj[0]=yj_initial
vjx[0]=vjx_initial
vjy[0]=vjy_initial

for i in range(npoints):
    #calulate the distance from asteroids to sun
    r1[i]=sqrt(x1[i]**2+y1[i]**2)
    r2[i]=sqrt(x2[i]**2+y2[i]**2)
    r3[i]=sqrt(x3[i]**2+y3[i]**2)
    r4[i]=sqrt(x4[i]**2+y4[i]**2)
    #the distance from Jupiter to sun
    rj[i]=sqrt(xj[i]**2+yj[i]**2)
    #the distance from asteroids to Jupiter
    r1j[i]=sqrt((x1[i]-xj[i])**2+(y1[i]-yj[i])**2)
    r2j[i]=sqrt((x2[i]-xj[i])**2+(y2[i]-yj[i])**2)
    r3j[i]=sqrt((x3[i]-xj[i])**2+(y3[i]-yj[i])**2)
    r4j[i]=sqrt((x4[i]-xj[i])**2+(y4[i]-yj[i])**2)
    #compute the x and y components for new velocity of earth
    v1x[i+1]=v1x[i]-4*pi**2*x1[i]*dt/r1[i]**3-4*pi**2*(Mj/Ms)*(x1[i]-xj[i])*dt/r1j[i]**3
    v1y[i+1]=v1y[i]-4*pi**2*y1[i]*dt/r1[i]**3-4*pi**2*(Mj/Ms)*(y1[i]-yj[i])*dt/r1j[i]**3
    #
    v2x[i+1]=v2x[i]-4*pi**2*x2[i]*dt/r2[i]**3-4*pi**2*(Mj/Ms)*(x2[i]-xj[i])*dt/r2j[i]**3
    v2y[i+1]=v2y[i]-4*pi**2*y2[i]*dt/r2[i]**3-4*pi**2*(Mj/Ms)*(y2[i]-yj[i])*dt/r2j[i]**3
    #
    v3x[i+1]=v3x[i]-4*pi**2*x3[i]*dt/r3[i]**3-4*pi**2*(Mj/Ms)*(x3[i]-xj[i])*dt/r3j[i]**3
    v3y[i+1]=v3y[i]-4*pi**2*y3[i]*dt/r3[i]**3-4*pi**2*(Mj/Ms)*(y3[i]-yj[i])*dt/r3j[i]**3
    #
    v4x[i+1]=v4x[i]-4*pi**2*x4[i]*dt/r4[i]**3-4*pi**2*(Mj/Ms)*(x4[i]-xj[i])*dt/r4j[i]**3
    v4y[i+1]=v4y[i]-4*pi**2*y4[i]*dt/r4[i]**3-4*pi**2*(Mj/Ms)*(y4[i]-yj[i])*dt/r4j[i]**3
    #compute the x and y components for new velocity of Jupiter
    vjx[i+1]=vjx[i]-4*pi**2*xj[i]*dt/rj[i]**3
    vjy[i+1]=vjy[i]-4*pi**2*yj[i]*dt/rj[i]**3
    #use Euler Cromer method to update the new positions of asteroids and Jupiter
    x1[i+1]=x1[i]+v1x[i+1]*dt
    y1[i+1]=y1[i]+v1y[i+1]*dt
    #
    x2[i+1]=x2[i]+v2x[i+1]*dt
    y2[i+1]=y2[i]+v2y[i+1]*dt
    #
    x3[i+1]=x3[i]+v3x[i+1]*dt
    y3[i+1]=y3[i]+v3y[i+1]*dt
    #
    x4[i+1]=x4[i]+v4x[i+1]*dt
    y4[i+1]=y4[i]+v4y[i+1]*dt
    #
    xj[i+1]=xj[i]+vjx[i+1]*dt
    yj[i+1]=yj[i]+vjy[i+1]*dt
    
fig=plt.figure(1)
plt.plot(x1,y1,'b.',x2,y2,'c.',x3,y3,'g.',x4,y4,'y.',xj,yj,'k')
plt.xlim(-5.5,5.5)
plt.ylim(-5.5,5.5)
plt.xlabel('x(AU)')
plt.ylabel('y(AU)')
plt.title('three body simulation-kirkwood gap&Jupiter')
fig.savefig('gap with ju')
fig=plt.figure(2)
plt.plot(x1,y1,'b.',x2,y2,'c.',x3,y3,'g.',x4,y4,'y.')
plt.xlim(-3.5,3.5)
plt.ylim(-3.5,3.5)
plt.xlabel('x(AU)')
plt.ylabel('y(AU)')
plt.title('three body simulation-kirkwood gap')
fig.savefig('gap')
plt.show()

