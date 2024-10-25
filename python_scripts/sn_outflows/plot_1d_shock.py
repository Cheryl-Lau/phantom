#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np 
import matplotlib.pyplot as plt 
import glob 
import os 


filename_confsn = 'semi_confined02/shock1d_withhii_chnl02_00200.dat'
filename_ffsn = 'free_field/shock1d_nocloud_bigenv_00020.dat'


# t1: shock1d_withhii_chnl02_00030.dat; shock1d_nocloud_bigenv_00003.dat
# t2: shock1d_withhii_chnl02_00050.dat; shock1d_nocloud_bigenv_00005.dat
# t3: shock1d_withhii_chnl02_00200.dat; shock1d_nocloud_bigenv_00020.dat
# t4: shock1d_withhii_chnl02_01400.dat; shock1d_nocloud_bigenv_00140.dat


#
# Semi-confined 
#

rho_scale = 1E-24
vel_scale = 1E5
u_scale = 1E14
thermpr_scale = 1E-9
rampr_scale = 1E-10

locx_cfshock = []
rho_cfshock = []
velx_cfshock = []
u_cfshock = []
thermpr_cfshock = []
rampr_cfshock = []

file = open(filename_confsn)
for iline,line in enumerate(file):
    line = line.split()
    if (iline == 1):
        time = line[0]
    elif (iline >= 3): 
        xloc = line[0]
        rho = line[1]
        velx = line[2]
        u = line[3]
        thermpr = line[4]
        rampr = line[5]
        locx_cfshock.append(float(xloc)/3.086E+18)
        rho_cfshock.append(float(rho)/rho_scale)
        velx_cfshock.append(float(velx)/vel_scale)
        u_cfshock.append(float(u)/u_scale)
        thermpr_cfshock.append(float(thermpr)/thermpr_scale)
        rampr_cfshock.append(float(rampr)/rampr_scale)

t0_cf = 3.34E+13 
time_Myr_cf = (float(time)-t0_cf)/(1E6*365*24*60*60)

#
# Free-field 
#

locx_ffshock = []
rho_ffshock = []
velx_ffshock = []
u_ffshock = []
thermpr_ffshock = []
rampr_ffshock = []

file = open(filename_ffsn)
for iline,line in enumerate(file):
    line = line.split()
    if (iline == 1):
        time = line[0]
    elif (iline >= 3): 
        xloc = line[0]
        rho = line[1]
        velx = line[2]
        u = line[3]
        thermpr = line[4]
        rampr = line[5]
        locx_ffshock.append(float(xloc)/3.086E+18)
        rho_ffshock.append(float(rho)/rho_scale)
        velx_ffshock.append(float(velx)/vel_scale)
        u_ffshock.append(float(u)/u_scale)
        thermpr_ffshock.append(float(thermpr)/thermpr_scale)
        rampr_ffshock.append(float(rampr)/rampr_scale)

time_Myr_ff = float(time)/(1E6*365*24*60*60)



fig1, (ax1,ax2) = plt.subplots(nrows=2, sharex=True, subplot_kw=dict(frameon=True), figsize=(8,10), dpi=300)
plt.subplots_adjust(hspace=.0)


ax1.set_ylim([8E-9,5E4])
ax1.set_yscale('log')
ax1.set_xlabel('x [pc]', fontsize=12)
ax1.tick_params(axis='x', labelsize=12)
ax1.tick_params(axis='y', labelsize=12)
ax1.text(0.05, 5e3, str(round(time_Myr_cf,3))+' Myr', fontsize=12)
ax1.text(20, 5e3, 'Semi-confined', fontsize=12)
ax1.plot(locx_cfshock,rho_cfshock,color='cyan',label='rho x 1E-24') # rho 
ax1.plot(locx_cfshock,velx_cfshock,color='lawngreen',label='vel x 1E+5')  # vel
ax1.plot(locx_cfshock,u_cfshock,color='gold',label='u x 1E+14')  # u
ax1.plot(locx_cfshock,thermpr_cfshock,color='magenta',label='P_therm x 1E-9')  # thermpr
ax1.plot(locx_cfshock,rampr_cfshock,color='darkorange',label='P_ram x 1E-9')  # rampr
ax1.legend()


ax2.set_ylim([8E-9,5E4])
ax2.set_yscale('log')
ax2.set_xlabel('x [pc]', fontsize=12)
ax2.tick_params(axis='x', labelsize=12)
ax2.tick_params(axis='y', labelsize=12)
ax2.text(0.05, 5e3, str(round(time_Myr_ff,3))+' Myr', fontsize=12)
ax2.text(20, 5e3, 'Free-field', fontsize=12)
ax2.plot(locx_cfshock,rho_ffshock,color='cyan',label='rho x 1E-24') # rho 
ax2.plot(locx_cfshock,velx_ffshock,color='lawngreen',label='vel x 1E+5')  # vel
ax2.plot(locx_cfshock,u_ffshock,color='gold',label='u x 1E+14')  # u
ax2.plot(locx_cfshock,thermpr_ffshock,color='magenta',label='P_therm x 1E-9')  # thermpr
ax2.plot(locx_cfshock,rampr_ffshock,color='darkorange',label='P_ram x 1E-9')  # rampr

fig1.savefig('shock_t3.png')

plt.show()

























