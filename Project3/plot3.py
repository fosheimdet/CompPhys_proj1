# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 10:16:53 2018

@author: foshe
"""

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import numpy as np
import matplotlib.pyplot as plt 



table1 = np.loadtxt("Mercury.txt",usecols= range(0), skiprows = 1)
table2 = np.loadtxt('Venus.txt',usecols=range(0),skiprows=1)
table3 = np.loadtxt('Earth.txt',usecols=range(0),skiprows=1)
table4 = np.loadtxt('Mars.txt',usecols=range(0),skiprows=1)
table5 = np.loadtxt('Jupiter.txt',usecols=range(0),skiprows=1)
table6 = np.loadtxt('Saturn.txt',usecols=range(0),skiprows=1)
table7 = np.loadtxt('Uranus.txt',usecols=range(0),skiprows=1)
table8 = np.loadtxt('Neptune.txt',usecols=range(0),skiprows=1)

table9 = np.loadtxt("precession.txt",usecols=range(0),skiprows=1 )
table10= np.loadtxt("MercuryRel.txt",usecols=range(0),skiprows=1 )


x1 = table1[:,0]; 
y1 = table1[:,1];

x2 = table2[:,0]; 
y2 = table2[:,1]; 

x3 = table3[:,0]; 
y3 = table3[:,1];

x4 = table4[:,0]; 
y4 = table4[:,1];

x5 = table5[:,0]; 
y5 = table5[:,1];

x6 = table6[:,0]; 
y6 = table6[:,1];

x7 = table7[:,0]; 
y7 = table7[:,1];

x8 = table8[:,0]; 
y8 = table8[:,1];

x9 = table9[:,0]
y9 = table9[:,1]

x10= table10[:,0]; 
y10= table10[:,0];



plt.figure(0)
plt.plot(x1,y1,label = "Mercury")
plt.plot(x2,y2,label = "Venus"); 
plt.plot(x3,y3, label = "Earth");
plt.plot(x4,y4, label = "Mars"); 
plt.plot(x5,y5, label = "Jupiter");
plt.plot(x6,y6, label = "Saturn"); 
plt.plot(x7,y7, label = "Uranus"); 
plt.plot(x8, y8, label = "Neptune" );  

plt.xlabel("x-position[AU]", size = 14)
plt.ylabel("y-position[AU]",size = 14) 
plt.title('Solar system (velocity verlet), N=1e6, t_max = 100yrs')
plt.legend()

plt.figure(1)
plt.plot(x9,y9)
plt.ylabel("Angle of perihelion", size = 14)
plt.xlabel("Number of orbits", size = 14)
plt.title('Precession of Mercurys perihelion, N=1e6, t_max = 100yrs')

plt.figure(2)
plt.plot(x10,x10)
plt.xlabel("x-position[AU]", size = 14)
plt.ylabel("y-position[AU]",size = 14)
plt.title('Orbit of Mercury with relativistic correction, N=1e6, t_max = 100yrs')
