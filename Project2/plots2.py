# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 11:40:48 2018

@author: foshe
"""
import numpy as np
import matplotlib.pyplot as plt 


#Loading all the data and storing them in tables
table0 = np.loadtxt('non_interacting.txt',usecols=range(0),skiprows=1)
table1 = np.loadtxt('omega0.01.txt',usecols=range(0),skiprows=1)
table2 = np.loadtxt('omega0.5.txt', usecols=range(0), skiprows=1)
table3 = np.loadtxt('omega1.txt', usecols=range(0), skiprows=1)
table4 = np.loadtxt('omega5.txt', usecols=range(0), skiprows=1)

table5 = np.loadtxt('bucklingbeam0.txt', usecols=range(0), skiprows = 1)
table6 = np.loadtxt('bucklingbeam1.txt', usecols=range(0), skiprows = 1)
table7 = np.loadtxt('bucklingbeam2.txt', usecols=range(0), skiprows = 1)
table8 = np.loadtxt('bucklingbeam3.txt', usecols=range(0), skiprows = 1)


rho = table0[:,0] #rho values for harmonic oscillator
u = table0[:,1]
v = table1[:,1]; 
w = table2[:,1]; 
x = table3[:,1]; 
y = table4[:,1]; 

Rho = table5[:,0] #rho values for buckling beam
a = table5[:,1];
b = table6[:,1];
c = table7[:,1];
d = table8[:,1];

#Finding probability functions for the different harmonic oscillator cases
u = np.abs(u)**2; 
v = np.abs(v)**2;
w = np.abs(w)**2;
x = np.abs(x)**2;
y = np.abs(y)**2;


plt.figure(0) 
plt.plot(rho,u,label="non-interacting")
plt.plot(rho,v,label="omega_r=0.01")
plt.legend()
plt.ylabel('|u(ρ)|^2', size=12)
plt.xlabel('ρ', size = 12)
plt.title('Harmonic oscillator with and without Coulomb potential')
#plt.xlim([0,5])
#plt.ylim([0,0.05])

plt.figure(1)
plt.plot(rho,v,label="omega_r=0.01")
plt.plot(rho,w,label="omega_r=0.5")
plt.plot(rho,x,label="omega_r=1")
plt.plot(rho,y,label="omega_r=5")
plt.legend()
plt.ylabel('|u(ρ)|^2', size = 14)
plt.xlabel('ρ', size = 14)
plt.title('Position probability for different potential strengths')
#plt.xlim([0,5])
#plt.ylim([0,0.1])

plt.figure(2) 
plt.plot(Rho,a,label="first eigenmode")
plt.plot(Rho,b,label="second eigenmode")
plt.plot(Rho,c,label="third eigenmode")
plt.plot(Rho,d,label="fourth eigenmode")
plt.legend()
plt.ylabel('u(ρ)', size=12)
plt.xlabel('ρ', size = 12)
plt.title('Eigenmodes of a buckling beam')
plt.xlim([0,1])
#plt.ylim([0,0.03])


