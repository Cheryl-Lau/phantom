#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 

plot_slice = True

data = np.loadtxt('xyzf_cminode_cropping.txt',skiprows=1)

radius = 5
x = data[:,0]
y = data[:,1]
z = data[:,2]
nH = data[:,3]

if (plot_slice == True): 
    x_slice = []
    y_slice = []
    nH_slice = []
    for i in range(len(x)):
        if z[i] > -0.1 and z[i] < 0.2:
            x_slice.append(x[i])
            y_slice.append(y[i])
            nH_slice.append(nH[i])

    x = x_slice 
    y = y_slice 
    nH = nH_slice

plt.xlim([-radius,radius])
plt.ylim([-radius,radius])

plt.scatter(x,y,s=1,c=nH,vmin=0.,vmax=1.)

plt.show() 







































