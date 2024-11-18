#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import glob 
import os 

pc = 3.086E+18
Myr = 1E6*365*24*60*60

# Sync to shock arrival time 
t0_shock_cf = 0.4E+14/Myr
t0_shock_ff = 0.03E+14/Myr


def plot_gas_detec(ax1,ax2,ax3,filepath,labelstr,colourstr,t0_shock):

    firstcall = True
    for filename in glob.glob(filepath):  

        time_cgs,x_cgs,y_cgs,z_cgs = np.loadtxt(filename,skiprows=1,max_rows=1,unpack=True)
        time_Myr = time_cgs/Myr - t0_shock

        rho_cgs,vx_cgs,u_cgs,thermpr_cgs,rampr_cgs = np.loadtxt(filename,skiprows=3,max_rows=3,unpack=True)

        if (firstcall==True):
            ax1.scatter(time_Myr,vx_cgs,s=3,color=colourstr,label=labelstr)
            ax2.scatter(time_Myr,thermpr_cgs,s=3,color=colourstr,label=labelstr)
            ax3.scatter(time_Myr,rampr_cgs,s=3,color=colourstr,label=labelstr)
            firstcall = False
        else:
            ax1.scatter(time_Myr,vx_cgs,s=3,color=colourstr)
            ax2.scatter(time_Myr,thermpr_cgs,s=3,color=colourstr)
            ax3.scatter(time_Myr,rampr_cgs,s=3,color=colourstr)

    return 

        
def plot_analyt_model(ax1,ax2,ax3,filepath,labelstr,colourstr):

    time_cgs,vx_cgs,thermpr_cgs,rampr_cgs = np.loadtxt(filepath,unpack=True)
    time_Myr = time_cgs/Myr 

    ax1.plot(time_Myr,vx_cgs,'--',color=colourstr,label=labelstr)
    ax2.plot(time_Myr,thermpr_cgs,'--',color=colourstr,label=labelstr)
    ax3.plot(time_Myr,rampr_cgs,'--',color=colourstr,label=labelstr)

    return 


fig,(ax1,ax2,ax3) = plt.subplots(nrows=3,sharex=False,subplot_kw=dict(frameon=True),figsize=(5,9),dpi=200)
#plt.subplots_adjust(hspace=.0)

plot_analyt_model(ax1,ax2,ax3,'semi_confined_model.txt','semi-confined model','salmon')
plot_analyt_model(ax1,ax2,ax3,'free_field_model.txt','free-field model','royalblue')

plot_gas_detec(ax1,ax2,ax3,'semi_confined01/gasflow_withhii_chnl01_*.dat','semi-confined','red',t0_shock_cf)
plot_gas_detec(ax1,ax2,ax3,'free_field/gasflow_nocloud_bigenv_*.dat','free-field','blue',t0_shock_ff)

ax1.set_xlabel('time [Myr]')
ax1.set_ylabel('velocity [cm $\mathrm{s^{-1}}$]')
ax1.set_yscale('log')
ax1.set_xlim([0,3.1])
ax1.set_ylim([1E5,6E7])
ax1.legend(loc='upper right',fontsize=8)

ax2.set_xlabel('time [Myr]')
ax2.set_ylabel('thermal pressure [g $\mathrm{cm^{-1}}$ $\mathrm{s^{-2}}$]')
ax2.set_yscale('log')
ax2.set_xlim([0,3.1])
ax2.set_ylim([2E-13,1E-8])
ax2.legend(loc='upper right',fontsize=8)

ax3.set_xlabel('time [Myr]')
ax3.set_ylabel('ram pressure [g $\mathrm{cm^{-1}}$ $\mathrm{s^{-2}}$]')
ax3.set_yscale('log')
ax3.set_xlim([0,3.1])
ax3.set_ylim([1E-19,5E-3])
ax3.legend(loc='upper right',fontsize=8)

fig.tight_layout(pad=1.0)

plt.savefig('gas_prop_detector.png')

plt.show()








