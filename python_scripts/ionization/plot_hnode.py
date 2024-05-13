# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 22:34:40 2023

@author: lausi
"""
import numpy as np 
import matplotlib.pyplot as plt 

x,y,z,h,m = np.loadtxt('xyzhm_cmi.txt',skiprows=1,unpack=True)

fig, ax = plt.subplots(figsize=(8, 8),dpi=100)
ax.set_aspect('equal', adjustable='box')

radius = 0.5
islice = np.where((z>-0.1)&(z<0.1))
x_slice = x[islice]
y_slice = y[islice]
h_slice = h[islice]

for i in range(len(x_slice)):
    circle = plt.Circle((x_slice[i],y_slice[i]),h_slice[i],alpha=0.2,color='blue',fill=False)
    ax.add_patch(circle)

ax.set_xlim([-radius,radius])
ax.set_ylim([-radius,radius])
ax.set_xlabel('x [pc]')
ax.set_ylabel('y [pc]')