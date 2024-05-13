#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import glob 
import os 


duration = 1.5E+13

#
# Confined case - Time evol of properties at target point 
#
time_evol_cav = []
rho_evol_cav = []
velx_evol_cav = []
thermpr_evol_cav = []
rampr_evol_cav = []

for filename in glob.glob('gasflow_cavitysn_*.dat'):  # si units 
    file = open(filename)
    for iline,line in enumerate(file):
        line = line.split()
        if (iline == 1):
            time = line[0]
            time_evol_cav.append(float(time))
        elif (iline == 3): 
            rho = line[0]
            velx = line[1]
            thermpr = line[2]
            rampr = line[3]
            rho_evol_cav.append(float(rho))
            velx_evol_cav.append(float(velx))
            thermpr_evol_cav.append(float(thermpr))
            rampr_evol_cav.append(float(rampr))

# Sync at point when sn shock arrives 
t0_cav = 3.411E+13
icrop = np.where(np.logical_and(np.array(time_evol_cav) > t0_cav, np.array(time_evol_cav) < t0_cav+duration))[0]

time_evol_cav = np.take(time_evol_cav,icrop) - t0_cav
rho_evol_cav = np.take(rho_evol_cav,icrop)
velx_evol_cav = np.take(velx_evol_cav,icrop)
thermpr_evol_cav = np.take(thermpr_evol_cav,icrop)
rampr_evol_cav = np.take(rampr_evol_cav,icrop)


#
# Free-field case - Time evol of properties at target point 
#
time_evol_ff = []
rho_evol_ff = []
velx_evol_ff = []
thermpr_evol_ff = []
rampr_evol_ff = []

for filename in glob.glob('gasflow_ffsn_*.dat'):  # si units 
    file = open(filename)
    for iline,line in enumerate(file):
        line = line.split()
        if (iline == 1):
            time = line[0]
            time_evol_ff.append(float(time))
        elif (iline == 3): 
            rho = line[0]
            velx = line[1]
            thermpr = line[2]
            rampr = line[3]
            rho_evol_ff.append(float(rho))
            velx_evol_ff.append(float(velx))
            thermpr_evol_ff.append(float(thermpr))
            rampr_evol_ff.append(float(rampr))

# Sync at point when sn shock arrives 
t0_ff = 6.2E+11 
icrop = np.where(np.logical_and(np.array(time_evol_ff) > t0_ff, np.array(time_evol_ff) < t0_ff+duration))[0]

time_evol_ff = np.take(time_evol_ff,icrop) - t0_ff
rho_evol_ff = np.take(rho_evol_ff,icrop)
velx_evol_ff = np.take(velx_evol_ff,icrop)
thermpr_evol_ff = np.take(thermpr_evol_ff,icrop)
rampr_evol_ff = np.take(rampr_evol_ff,icrop)



'''
fig1 = plt.figure(figsize=[7,5])
ax1 = fig1.add_subplot(111)
ax1.scatter(time_evol_cav,rho_evol_cav,s=1,color='red',label='partially-confined')
ax1.scatter(time_evol_ff,rho_evol_ff,s=1,color='blue')
ax1.set_xlabel('time [s]')
ax1.set_ylabel('rho [kg m^-3]')
ax1.set_yscale('log')
ax1.legend()
'''

fig2 = plt.figure(figsize=[7,5])
ax2 = fig2.add_subplot(111)
ax2.scatter(time_evol_cav,velx_evol_cav,s=1,color='red',label='partially-confined')
ax2.scatter(time_evol_ff,velx_evol_ff,s=1,color='blue',label='free-field')
ax2.set_xlabel('time [s]')
ax2.set_ylabel('v_x [m s^-1]')
ax2.set_yscale('log')
ax2.set_ylim([5E3,6E5])
ax2.legend()
fig2.savefig('sim_vx.png',dpi=200)

fig3 = plt.figure(figsize=[7,5])
ax3 = fig3.add_subplot(111)
ax3.scatter(time_evol_cav,thermpr_evol_cav,s=1,color='red',label='partially-confined')
ax3.scatter(time_evol_ff,thermpr_evol_ff,s=1,color='blue',label='free-field')
ax3.set_xlabel('time [s]')
ax3.set_ylabel('therm pr [kg m^-1 s^-2]')
ax3.set_yscale('log')
ax3.set_ylim([2E-13,6E-6])
ax3.legend()
fig3.savefig('sim_thermp.png',dpi=200)

fig4 = plt.figure(figsize=[7,5])
ax4 = fig4.add_subplot(111)
ax4.scatter(time_evol_cav,rampr_evol_cav,s=1,color='red',label='partially-confined')
ax4.scatter(time_evol_ff,rampr_evol_ff,s=1,color='blue',label='free-field')
ax4.set_xlabel('time [s]')
ax4.set_ylabel('ram pr [kg m^-1 s^-2]')
ax4.set_yscale('log')
ax4.set_ylim([3E-15,5E-5])
ax4.legend()
fig4.savefig('sim_ramp.png',dpi=200)

plt.show()





















