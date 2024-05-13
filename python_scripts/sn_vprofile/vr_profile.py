# -*- coding: utf-8 -*-


import numpy as np 
import matplotlib.pyplot as plt 


r_dist = 0.05
r_sn = 0.2
r = 2*r_dist/r_sn


def vrfunc(a,r): 
    if r >= 0 and r < 1: 
        vr = a*(-3*r+9/4*r**2) 
    elif r >= 1 and r < 2: 
        vr = -3/4*a*(2-r)**2
    else:
        vr = 0
    return -vr

a = 1

r_array = np.linspace(0,2,100) 
vr_array = [vrfunc(a,r) for r in r_array]

plt.rcParams['font.size'] = '12'

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
ax.plot(r_array,vr_array,color='black')
ax.set_xlabel('r',fontsize=12, rotation=0, labelpad=10)
ax.set_ylabel(r'$ \mathrm{ v_r } $',fontsize=12, rotation=90, labelpad=10)
ax.tick_params(axis='x', direction='in')
ax.tick_params(axis='y', direction='in')
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
               bottom=True, top=True, left=True, right=True)

ax.tick_params(axis='both', which='major', pad=10)

plt.show()




































