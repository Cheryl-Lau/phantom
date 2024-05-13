# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 20:54:53 2023

@author: lausi
"""

import numpy as np 
import matplotlib.pyplot as plt 

n,i,x,y,z,h,m,f = np.loadtxt('nixyzhmf_twosrc_lowresfront.txt',unpack=True,skiprows=3)

radius = 0.7
plt.figure(figsize=[6,5],dpi=200)
plt.xlim([-radius,radius])
plt.ylim([-radius,radius])
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

islice = np.where((z>-0.07)&(z<0.07))
plt.scatter(x[islice],y[islice],s=0.5,c=f[islice],cmap='viridis')
cb = plt.colorbar(orientation='vertical', shrink=0.85)
cb.set_label('Neutral fraction') 
plt.xlabel('x [pc]')
plt.ylabel('y [pc]')
plt.show()

plt.savefig('node_twosource_nH.png')