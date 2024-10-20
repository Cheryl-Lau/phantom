#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import glob 
import os 

pc = 3.086E+18
Myr = 1E6*365*24*60*60

def extract_shell_rad(ax):

    t0_shock = 1.075  # Myr 
    firstcall = True

    for filename in glob.glob('semi_confined01/shell_withhii_chnl01_*.dat'):  # cgs units 

        time_cgs = np.loadtxt(filename,skiprows=1,max_rows=1)
        time_Myr = time_cgs/Myr - t0_shock

        x_cgs,rho_cgs,u_cgs = np.loadtxt(filename,skiprows=3,unpack=True)
        x_pc = x_cgs/pc

        # Extract the large jumps off cliff in log(u) plot 
        grads = np.gradient(np.log10(u_cgs))
        ineg = np.where(grads < 0)[0]
        largegrad = np.percentile(np.abs(grads[ineg]),70)
        ilargeneg = np.where(grads < -1*largegrad)
        r_shell_pc_u = x_pc[ilargeneg]
        time_Myr_arr_u = np.full(len(r_shell_pc_u),time_Myr)

        # Extract the large jumps off cliff in log(rho) plot 
        grads = np.gradient(np.log10(rho_cgs))
        ineg = np.where(grads < 0)[0]
        largegrad = np.percentile(np.abs(grads[ineg]),80)
        ilargeneg = np.where(grads < -1*largegrad)
        r_shell_pc_rho = x_pc[ilargeneg]
        time_Myr_arr_rho = np.full(len(r_shell_pc_rho),time_Myr)

        # Plot both 
        if (firstcall==True):
            ax.scatter(time_Myr_arr_u,np.abs(r_shell_pc_u),s=1,color='lightskyblue',label='large '+r'$-\nabla u$')
            ax.scatter(time_Myr_arr_rho,np.abs(r_shell_pc_rho),s=1,color='salmon',label='large '+r'$-\nabla \rho$')
            firstcall = False
        else:
            ax.scatter(time_Myr_arr_u,np.abs(r_shell_pc_u),s=1,color='lightskyblue')
            ax.scatter(time_Myr_arr_rho,np.abs(r_shell_pc_rho),s=1,color='salmon')


def plot_analyt_model(ax):

    time_cgs,rad_shell_cgs,vel_shell_cgs = np.loadtxt('shell_evol_model.txt',unpack=True)
    time_Myr = time_cgs/Myr 
    rad_shell_pc = rad_shell_cgs/pc 

    ax.plot(time_Myr,rad_shell_pc,color='black',label='analytical model')


fig = plt.figure(figsize=[7,5],dpi=200)
ax = fig.add_subplot(111)    

plot_analyt_model(ax)
extract_shell_rad(ax)

ax.set_xlabel('time [Myr]')
ax.set_ylabel('cavity radius [pc]')
ax.set_yscale('log')
ax.set_ylim([2,280])
ax.legend()

plt.savefig('shell_radius_evol.png')
plt.show()












