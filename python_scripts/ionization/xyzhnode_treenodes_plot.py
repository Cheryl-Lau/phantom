#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 

plot_slice = True
colour_mass = False  # else nH

time_cgs = np.loadtxt('nixyzhmf_04999.txt',skiprows=1,max_rows=1)
Myr = 1e6*365*24*60*60 

n,i,x,y,z,h,m,nH = np.loadtxt('nixyzhmf_04999.txt',skiprows=3,unpack=True)
imfilename = 'nixyzhmf_04999.png'

radius = 4.5
centre = [2.181616076E+00,2.883273141E+00,-2.588239455E+00]  # sink 139
zmin = -2.5
zmax = -1.0

print(np.sum(m))

max_m = max(m)
m_scaled = m/max_m

if (plot_slice == True): 

    islice = np.where(np.logical_and(z > zmin,z < zmax))
    x_slice = x[islice]
    y_slice = y[islice]
    m_slice = m[islice]
    h_slice = h[islice]
    nH_slice = nH[islice]

    max_m = max(m_slice)
    m_scaled = m_slice/max_m

    x = x_slice 
    y = y_slice 
    nH = nH_slice


fig, ax1 = plt.subplots(ncols=1, sharey=True, subplot_kw=dict(frameon=True), figsize=(6,5))

if (colour_mass == True):
    m = ax1.scatter(x,y,s=1,c=m_scaled,alpha=0.6)
else:
    m = ax1.scatter(x,y,s=1,c=nH,alpha=0.6)
    
plt.gca().set_aspect('equal')

ax1.set_xlim([-radius+centre[0],radius+centre[0]])
ax1.set_ylim([-radius+centre[1],radius+centre[1]])
ax1.set_xlabel('x [pc]')
ax1.set_ylabel('y [pc]')
ax1.text(-radius*0.9+centre[0],radius*0.85+centre[1],'time = '+str(round(time_cgs/Myr,2))+' Myr')
ax1.text(radius*0.2+centre[0],radius*0.85+centre[1],str(zmin)+' pc < z < '+str(zmax)+' pc')

plt.colorbar(m,label='neutral fraction', orientation='vertical') 


fig.tight_layout(pad=0.5)


plt.savefig(imfilename.strip(),dpi=200)
plt.show() 







































