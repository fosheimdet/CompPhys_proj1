# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 18:19:19 2018

@author: foshe
"""

import numpy as np
import matplotlib.pyplot as plt 

#f=open("data.txt","r")
#if f.mode == "r": 
#    contents = f.read()
#    print(contents)
#f.close()

table = np.loadtxt('gaussian_elimination1000x1000.txt',usecols=range(0),skiprows=1)

x = table[:,0]
u = table[:,1]
v = table[:,2]
print(v)

plt.plot(x,u,label="Analytical solution")
plt.plot(x,v,'--',label="Numerical solution")
plt.legend()
plt.ylabel('u,v')
plt.xlabel('x')
plt.title('n=1000, h=9.99*10^-4')

plt.ylim(0,0.8)
plt.xlim(0,1)

#for i in range (1,20):
#    print(table(i,1))

#with open("data.txt", "r") as f: 
#    
#    for i in range (1,10):
#        f_contents = f.read(8)
#        f_blank = f.read(1)
#        f_contents2=f.read(8)
#        f_blank = f.read(1)
##        print(f_contents,"\n")
#        data = float(f_contents)
#        
#        print(data)
    

    
 

