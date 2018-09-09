# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 15:19:49 2018

@author: foshe
"""


import matplotlib.pyplot as plt 
logh = [-1,-2,-3,-4,-5,-6,-7, -7.2]
e = [-1.1797, -3.08804, -5.08005, -7.07936, - 9.0049, -6.77137, -13.007, -7.25184]


plt.plot(logh, e)
plt.show()
plt.ylabel('max(log_10(|v_i-u_i)/u_i|))')
plt.xlabel('log_10(h)')
plt.title('Relative error vs steplength')
