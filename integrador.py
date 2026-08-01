#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 14:25:29 2017

Integra un sistema 2D y dibuja retrato de fases

@author: juan
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#%%
def f(z,t):
    #z = [x, y]
    x = z[0]
    y = z[1]
    #return([dx/dt,dy/dt])
    return([np.sin(y) ,x-x**3])

tmax = 3.	# tiempo de integración hacia adelante
dt = 1./100000
t = np.linspace(0, tmax, int(tmax/dt))
#Figura
fig = plt.figure()
ax = fig.add_subplot(111)
#Condiciones iniciales
Xi = np.linspace(-1.5,1.5,3)
Yi = np.linspace(-1,5,8)

# Integro para cada condicion inicial
for i in Xi:
    for j in Yi:
        xi = [i,j]
        sol = odeint(f, xi, t)	# Integro
        #Concateno las integraciones
        x = sol[:,0]
        y = sol[:,1]
        ax.plot(x,y)

#Flechas de flujo
X = np.linspace(-10,10,30)
Y = np.linspace(-10,10,30)
XX,YY = np.meshgrid(X,Y)
DX,DY = f([XX,YY],t)
M = (np.hypot(DX, DY))
# Normalizo longitud
M[ M == 0] = 1.
DX /= M 
DY /= M 
ax.quiver(XX,YY,DX,DY, pivot = 'mid')

ax.set_xlim([-3,3])
ax.set_ylim([-10,10])