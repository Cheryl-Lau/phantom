#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import glob 
import os 

pc = 3.086E+18
Myr = 1E6*365*24*60*60

# Sync to shock arrival time 
t0_shock_scf = 0.091
t0_shock_cf = 0.172 
t0_shock_ff = 0.023


def plot_gas_detec(ax1,ax2,ax3,filepath,labelstr,colourstr,t0_shock):

    firstcall = True
    for filename in glob.glob(filepath):  

        time_cgs,x_cgs,y_cgs,z_cgs = np.loadtxt(filename,skiprows=1,max_rows=1,unpack=True)
        time_Myr = time_cgs/Myr - t0_shock

        rho_cgs,vx_cgs,u_cgs,thermpr_cgs,rampr_cgs = np.loadtxt(filename,skiprows=3,max_rows=3,unpack=True)

        if (firstcall==True):
            ax1.scatter(time_Myr,vx_cgs,s=2,color=colourstr,label=labelstr)
            ax2.scatter(time_Myr,thermpr_cgs,s=2,color=colourstr,label=labelstr)
            ax3.scatter(time_Myr,rampr_cgs,s=2,color=colourstr,label=labelstr)
            firstcall = False
        else:
            ax1.scatter(time_Myr,vx_cgs,s=2,color=colourstr)
            ax2.scatter(time_Myr,thermpr_cgs,s=2,color=colourstr)
            ax3.scatter(time_Myr,rampr_cgs,s=2,color=colourstr)

    return 

        
def plot_analyt_model(ax1,ax2,ax3,filepath,labelstr,colourstr):

    time_cgs,vx_cgs,thermpr_cgs,rampr_cgs = np.loadtxt(filepath,unpack=True)
    time_Myr = time_cgs/Myr 

    ax1.plot(time_Myr,vx_cgs,'--',color=colourstr,label=labelstr)
    ax2.plot(time_Myr,thermpr_cgs,'--',color=colourstr,label=labelstr)
    ax3.plot(time_Myr,rampr_cgs,'--',color=colourstr,label=labelstr)

    return 


fig,(ax1,ax2,ax3) = plt.subplots(nrows=3,sharex=False,subplot_kw=dict(frameon=True),figsize=(5,9),dpi=200)

plot_analyt_model(ax1,ax2,ax3,'semi_confined_model.txt','semi-confined model','salmon')
plot_analyt_model(ax1,ax2,ax3,'free_field_model.txt','free-field model','lightblue')
plot_analyt_model(ax1,ax2,ax3,'confined_model.txt','confined model','lightgreen')

plot_gas_detec(ax1,ax2,ax3,'semi_confined/vary_chnlsize/chnl005/gasflow_semiconf_1chnl005_*.dat','semi-confined','red',t0_shock_scf)
plot_gas_detec(ax1,ax2,ax3,'free_field/gasflow_freefield_*.dat','free-field','blue',t0_shock_ff)
plot_gas_detec(ax1,ax2,ax3,'confined/gasflow_confined_*.dat','confined','green',t0_shock_cf)

ax1.set_xlabel('time [Myr]')
ax1.set_ylabel('velocity [cm $\mathrm{s^{-1}}$]')
ax1.set_yscale('log')
ax1.set_xlim([0,0.38])
ax1.set_ylim([2E5,8E8])
ax1.legend(loc='upper right',fontsize=8)

ax2.set_xlabel('time [Myr]')
ax2.set_ylabel('thermal pressure [g $\mathrm{cm^{-1}}$ $\mathrm{s^{-2}}$]')
ax2.set_yscale('log')
ax2.set_xlim([0,0.38])
ax2.set_ylim([6E-13,2E-4])
ax2.legend(loc='upper right',fontsize=8)

ax3.set_xlabel('time [Myr]')
ax3.set_ylabel('ram pressure [g $\mathrm{cm^{-1}}$ $\mathrm{s^{-2}}$]')
ax3.set_yscale('log')
ax3.set_xlim([0,0.38])
ax3.set_ylim([1E-19,9E2])
ax3.legend(loc='upper right',fontsize=8)

fig.tight_layout(pad=1.0)

plt.savefig('gas_prop_detector_adiabatic.pdf',format='pdf')
plt.savefig('gas_prop_detector_adiabatic.png')


plt.show()





