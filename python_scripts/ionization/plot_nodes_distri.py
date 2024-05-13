# -*- coding: utf-8 -*-
"""
Plots the distribution of nodes for a given file
colour for mass or nH
"""
import numpy as np
import matplotlib.pyplot as plt 
import os


plot_tree = True  # else parts
plot_mass = True  # else nH
plot_h = True

pwd = os.getcwd()
pardir = os.path.abspath(os.path.join(pwd, os.pardir))

if plot_tree: 
    filename = pardir+'\\tree_snapshots\\nixyzhmf_00001.txt'
    nixyzhmf = np.loadtxt(filename,skiprows=3)
    x = nixyzhmf[:,2]
    y = nixyzhmf[:,3]
    z = nixyzhmf[:,4]
    h = nixyzhmf[:,5]
    m = nixyzhmf[:,6]
    nH = nixyzhmf[:,7]
    
    label = "nodes_"
    for c in filename:
        if c.isdigit():
            label = label + c
else:
    filename = pardir+'\\parts_snapshots\\xyzhmf_00090.txt'
    xyzhmf = np.loadtxt(filename,skiprows=3)
    x = xyzhmf[:,0]
    y = xyzhmf[:,1]
    z = xyzhmf[:,2]
    h = xyzhmf[:,3]
    m = xyzhmf[:,4]
    nH = xyzhmf[:,5]
    
    label = "parts_"
    for c in filename:
        if c.isdigit():
            label = label + c
    
islice = np.where(np.logical_and(z > -0.02,z < 0.02))
x_slice = x[islice]
y_slice = y[islice]
m_slice = m[islice]
h_slice = h[islice]
nH_slice = nH[islice]

radius = 0.9
plt.figure(figsize=[10,10],dpi=200)
plt.xlim([-radius,radius])
plt.ylim([-radius,radius])
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

if plot_mass == True:
    plt.scatter(x_slice,y_slice,s=1,c=m_slice,cmap='plasma',vmin=0.,vmax=0.035)
    cb = plt.colorbar(orientation='vertical', shrink=0.75)
    cb.set_label('Node / Particle mass')
    label = label + '_mass'
else:
    plt.scatter(x_slice,y_slice,s=1,c=nH_slice,cmap='viridis')
    cb = plt.colorbar(orientation='vertical', shrink=0.75)
    cb.set_label('Neutral fraction')
    label = label + '_nH'
    
plt.xlabel('x [pc]')
plt.ylabel('y [pc]')
    
plt.show()
plt.savefig(label+'.png')

if plot_h == True: 
    fig, ax = plt.subplots(figsize=(8, 8),dpi=200)
    ax.set_aspect('equal', adjustable='box')

    for i in range(len(x_slice)):
        circle = plt.Circle((x_slice[i],y_slice[i]),h_slice[i],alpha=0.2,color='blue',fill=False)
        ax.add_patch(circle)

    ax.set_xlim([-radius,radius])
    ax.set_ylim([-radius,radius])
    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    
    fname = 'hnode'
    if plot_tree == True:
        fname = fname + '_nodes'
    else:
        fname = fname + '_parts'
    fig.savefig(fname+'.png')

