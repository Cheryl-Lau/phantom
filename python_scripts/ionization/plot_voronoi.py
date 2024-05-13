# -*- coding: utf-8 -*-
"""
Draws voronoi grids around particles or nodes 
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d



data = np.loadtxt('nixyzhmf_cminode_twosource.txt',skiprows=1)
x = data[:,2]
y = data[:,3]
z = data[:,4]

points = []
for i in range(len(x)):
    if z[i] > -0.03 and z[i] < 0.0:
        if x[i] > -0.4 and x[i] < 0.4:
            if y[i] > -0.4 and y[i] < 0.6:
                points.append([x[i],y[i]])
points = np.array(points)
print(np.shape(points))

vor = Voronoi(points)

fig = voronoi_plot_2d(vor)

fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange',
                      line_width=2, line_alpha=0.6, point_size=2)
plt.show()

plt.figure()
plt.scatter(points[:,0],points[:,1],s=5,color='black')
plt.show()

