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

        # Plot both 
        if (firstcall==True):
            time_Myr_arr,r_shell_pc = get_large_falls_norm(u_cgs,70,time_Myr,x_pc)
            ax.scatter(time_Myr_arr,np.abs(r_shell_pc),s=1,color='paleturquoise',label='$70^\mathrm{th}$ percentile of '+r'$-\nabla u$')
            time_Myr_arr,r_shell_pc = get_large_falls_norm(u_cgs,80,time_Myr,x_pc)
            ax.scatter(time_Myr_arr,np.abs(r_shell_pc),s=1,color='lightskyblue',label='$80^\mathrm{th}$ percentile of '+r'$-\nabla u$')
            time_Myr_arr,r_shell_pc = get_large_falls_norm(u_cgs,90,time_Myr,x_pc)
            ax.scatter(time_Myr_arr,np.abs(r_shell_pc),s=1,color='cornflowerblue',label='$90^\mathrm{th}$ percentile of '+r'$-\nabla u$')

            time_Myr_arr,r_shell_pc = get_large_falls(rho_cgs,70,time_Myr,x_pc)
            ax.scatter(time_Myr_arr,np.abs(r_shell_pc),s=1,color='mistyrose',label='$70^\mathrm{th}$ percentile of '+r'$-\nabla \rho$')
            time_Myr_arr,r_shell_pc = get_large_falls(rho_cgs,80,time_Myr,x_pc)
            ax.scatter(time_Myr_arr,np.abs(r_shell_pc),s=1,color='lightsalmon',label='$80^\mathrm{th}$ percentile of '+r'$-\nabla \rho$')
            time_Myr_arr,r_shell_pc = get_large_falls(rho_cgs,90,time_Myr,x_pc)
            ax.scatter(time_Myr_arr,np.abs(r_shell_pc),s=1,color='coral',label='$90^\mathrm{th}$ percentile of '+r'$-\nabla \rho$')
            firstcall = False

        else:
            time_Myr_arr,r_shell_pc = get_large_falls_norm(u_cgs,70,time_Myr,x_pc)
            ax.scatter(time_Myr_arr,np.abs(r_shell_pc),s=1,color='paleturquoise')
            time_Myr_arr,r_shell_pc = get_large_falls_norm(u_cgs,80,time_Myr,x_pc)
            ax.scatter(time_Myr_arr,np.abs(r_shell_pc),s=1,color='lightskyblue')
            time_Myr_arr,r_shell_pc = get_large_falls_norm(u_cgs,90,time_Myr,x_pc)
            ax.scatter(time_Myr_arr,np.abs(r_shell_pc),s=1,color='cornflowerblue')

            time_Myr_arr,r_shell_pc = get_large_falls(rho_cgs,70,time_Myr,x_pc)
            ax.scatter(time_Myr_arr,np.abs(r_shell_pc),s=1,color='mistyrose')
            time_Myr_arr,r_shell_pc = get_large_falls(rho_cgs,80,time_Myr,x_pc)
            ax.scatter(time_Myr_arr,np.abs(r_shell_pc),s=1,color='lightsalmon')
            time_Myr_arr,r_shell_pc = get_large_falls(rho_cgs,90,time_Myr,x_pc)
            ax.scatter(time_Myr_arr,np.abs(r_shell_pc),s=1,color='coral')



def get_large_falls(physprop_cgs,percentile,time_Myr,x_pc):

    grads = np.gradient(np.log10(physprop_cgs))#/np.log10(physprop_cgs)
    ineg = np.where(grads < 0)[0]
    largegrad = np.percentile(np.abs(grads[ineg]),percentile)
    ilargeneg = np.where(grads < -1*largegrad)
    r_shell_pc = x_pc[ilargeneg]
    time_Myr_arr = np.full(len(r_shell_pc),time_Myr)

    return time_Myr_arr,r_shell_pc


def get_large_falls_norm(physprop_cgs,percentile,time_Myr,x_pc):

    grads = np.gradient(np.log10(physprop_cgs))/np.log10(physprop_cgs)
    ineg = np.where(grads < 0)[0]
    largegrad = np.percentile(np.abs(grads[ineg]),percentile)
    ilargeneg = np.where(grads < -1*largegrad)
    r_shell_pc = x_pc[ilargeneg]
    time_Myr_arr = np.full(len(r_shell_pc),time_Myr)

    return time_Myr_arr,r_shell_pc


def plot_analyt_model(ax):

    time_cgs,rad_shell_cgs,vel_shell_cgs = np.loadtxt('shell_evol_model.txt',unpack=True)
    time_Myr = time_cgs/Myr 
    rad_shell_pc = rad_shell_cgs/pc 

    ax.plot(time_Myr,rad_shell_pc,color='black',label='analytical model')


fig = plt.figure(figsize=[5.5,4],dpi=200)
ax = fig.add_subplot(111)    

plot_analyt_model(ax)
extract_shell_rad(ax)

ax.set_xlabel('time [Myr]')
ax.set_ylabel('cavity radius [pc]')
ax.set_yscale('log')
ax.set_ylim([0.9,150])
ax.legend(loc='lower right',fontsize=8)

fig.tight_layout(pad=0.5)

plt.savefig('shell_radius_evol.png')
plt.show()












