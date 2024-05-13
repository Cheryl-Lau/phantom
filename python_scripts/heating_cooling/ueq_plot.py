#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot ueqs wrt rho and gamma 
"""

import numpy as np 
import matplotlib.pyplot as plt 



data = np.loadtxt('rho_gamma_ueqs.txt',skiprows=1)

rho = data[:,0]
gamma = data[:,1]
ueq1 = data[:,3]
ueq2 = data[:,4]
ueq3 = data[:,5] 

# stores all ind points 
rho_all = []
gamma_all = []
ueq_all = []

for i in range(len(rho)):
    if (ueq1[i] != 0.):
        rho_all.append(rho[i])
        gamma_all.append(gamma[i])
        ueq_all.append(ueq1[i])
        
for i in range(len(rho)):
    if (ueq2[i] != 0.):
        rho_all.append(rho[i])
        gamma_all.append(gamma[i])
        ueq_all.append(ueq2[i])

for i in range(len(rho)):
    if (ueq3[i] != 0.):
        rho_all.append(rho[i])
        gamma_all.append(gamma[i])
        ueq_all.append(ueq3[i])
        
# 3D plot of all points 
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(np.log10(rho_all), np.log10(gamma_all), np.log10(ueq_all))


ax.set_xlabel('rho')
ax.set_ylabel('gamma')
ax.set_zlabel('u_eq')

plt.show()

