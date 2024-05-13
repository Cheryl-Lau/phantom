#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 

plot_slice = True
colour_mass = False  # else nH

data = np.loadtxt('nixyzhmf_cminode.txt',skiprows=1)

radius = 5
x = data[:,2]
y = data[:,3]
z = data[:,4]
m = data[:,6]
nH = data[:,7]

print(np.sum(m))

max_m = max(m)
m_scaled = m/max_m

if (plot_slice == True): 
    x_slice = []
    y_slice = []
    m_slice = []
    nH_slice = []
    for i in range(len(x)):
        if z[i] > -0.1 and z[i] < 0.2:
            x_slice.append(x[i])
            y_slice.append(y[i])
            m_slice.append(m[i])
            nH_slice.append(nH[i])
    max_m = max(m_slice)
    m_scaled = m_slice/max_m

    x = x_slice 
    y = y_slice 
    nH = nH_slice

plt.xlim([-radius,radius])
plt.ylim([-radius,radius])
if (colour_mass == True):
    plt.scatter(x,y,s=1,c=m_scaled)
else:
    plt.scatter(x,y,s=1,c=nH)
    
plt.show() 







































