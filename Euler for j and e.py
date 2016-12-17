# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 20:25:08 2016

@author: xh810
"""

import numpy as np
from math import pi,sqrt
import matplotlib.pyplot as plt

npoints=1000000
dt=0.001
Ms=2.0e30 #mass of sun in kg
Me=6.0e24 #mass of earth in kg
Mj=1.9e27 #mass of jupiter in kg

xe_initial=1 #initial position of earth in AU
ye_initial=0
vex_initial=0 #initial velocity of earth in AU/yr
vey_initial=2*pi

xj_initial=-5.2 #initial position of Jupiter in AU; assume at opposition initially
yj_initial=0
vjx_initial=0  #iniital velocity of Jupiter
vjy_initial=-2.7549 #2*pi*5.2 AU/11.85 years=2.7549 AU/yr

#create array to store position and velocity of earth
xe=np.zeros(npoints+1,float)
ye=np.zeros(npoints+1,float)
vex=np.zeros(npoints+1,float)
vey=np.zeros(npoints+1,float)

#create array to store position and velocity of Jupoter
xj=np.zeros(npoints+1,float)
yj=np.zeros(npoints+1,float)
vjx=np.zeros(npoints+1,float)
vjy=np.zeros(npoints+1,float)

re=np.zeros(npoints+1,float)
rj=np.zeros(npoints+1,float)
rej=np.zeros(npoints+1,float)

#initial position and velocity of earth
xe[0]=xe_initial
ye[0]=ye_initial
vex[0]=vex_initial
vey[0]=vey_initial

#initial position and velocity of Jupiter
xj[0]=xj_initial
yj[0]=yj_initial
vjx[0]=vjx_initial
vjy[0]=vjy_initial

for i in range(npoints):
    #calulate the distance from earth to sun
    re[i]=sqrt(xe[i]**2+ye[i]**2)
    #the distance from Jupiter to sun
    rj[i]=sqrt(xj[i]**2+yj[i]**2)
    #the distance from earth to Jupiter
    rej[i]=sqrt((xe[i]-xj[i])**2+(ye[i]-yj[i])**2)
    #compute the x and y components for new velocity of earth
    vex[i+1]=vex[i]-4*pi**2*xe[i]*dt/re[i]**3-4*pi**2*(Mj/Ms)*(xe[i]-xj[i])*dt/rej[i]**3
    vey[i+1]=vey[i]-4*pi**2*ye[i]*dt/re[i]**3-4*pi**2*(Mj/Ms)*(ye[i]-yj[i])*dt/rej[i]**3
    #compute the x and y components for new velocity of Jupiter
    vjx[i+1]=vjx[i]-4*pi**2*xj[i]*dt/rj[i]**3-4*pi**2*(Me/Ms)*(xj[i]-xe[i])*dt/rej[i]**3
    vjy[i+1]=vjy[i]-4*pi**2*yj[i]*dt/rj[i]**3-4*pi**2*(Me/Ms)*(yj[i]-ye[i])*dt/rej[i]**3
    #use Euler Cromer method to update the new positions of earth and Jupiter
    xe[i+1]=xe[i]+vex[i+1]*dt
    ye[i+1]=ye[i]+vey[i+1]*dt
    xj[i+1]=xj[i]+vjx[i+1]*dt
    yj[i+1]=yj[i]+vjy[i+1]*dt
    
fig=plt.figure()
plt.plot(xe,ye,'r',xj,yj,'k')
plt.xlim(-7,7)
plt.ylim(-7,7)
plt.xlabel('x(AU)')
plt.ylabel('y(AU)')
plt.title('three body simulation-Jupiter&Earth')
fig.savefig('jupiter and earth')
plt.show()

