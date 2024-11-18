#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


pc = 3.086E+18
Myr = 1E6*365*24*60*60



def get_large_falls(ax,physprop_cgs,percentile,time_Myr,x_pc):

    grads = np.gradient(np.log10(physprop_cgs))
    ineg = np.where(grads < 0)[0]
    largegrad = np.percentile(np.abs(grads[ineg]),percentile)
    ilargeneg = np.where(grads < -1*largegrad)
    r_shell_pc = x_pc[ilargeneg]
    time_Myr_arr = np.full(len(r_shell_pc),time_Myr)


    return time_Myr_arr,r_shell_pc



def main():

    filename = 'semi_confined01/shell_withhii_chnl01_01400.dat'
    t0_shock = 1.075  # Myr 

    time_cgs = np.loadtxt(filename,skiprows=1,max_rows=1)
    time_Myr = time_cgs/Myr - t0_shock

    x_cgs,rho_cgs,u_cgs = np.loadtxt(filename,skiprows=3,unpack=True)
    x_pc = x_cgs/pc


    fig = plt.figure(figsize=(10,10),dpi=200)
    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(323)
    ax3 = fig.add_subplot(325)
    ax4 = fig.add_subplot(322)
    ax5 = fig.add_subplot(324)
    ax6 = fig.add_subplot(326)

    ax1.plot(x_pc,rho_cgs)
    ax2.plot(x_pc,-1*np.gradient(rho_cgs))
    ax3.plot(x_pc,-1*np.gradient(rho_cgs)/rho_cgs)

    ax4.plot(x_pc,u_cgs)
    ax5.plot(x_pc,-1*np.gradient(u_cgs))
    ax6.plot(x_pc,-1*np.gradient(u_cgs)/u_cgs)


    ax1.set_ylabel(r'$\rho$')
    ax1.set_yscale('log')
    ax2.set_ylabel(r'$-\nabla \rho$')
    ax2.set_yscale('log')
    ax3.set_ylabel(r'$-\nabla \rho / \rho$')
    ax3.set_yscale('log')
    ax3.set_xlabel('x [pc]')

    ax4.set_ylabel(r'$u$')
    ax4.set_yscale('log')
    ax5.set_ylabel(r'$-\nabla u$')
    ax5.set_yscale('log')
    ax6.set_ylabel(r'$-\nabla u / u$')
    ax6.set_yscale('log')
    ax6.set_xlabel('x [pc]')

    ax1.text(-70,1e-22,str(round(time_Myr,2))+' Myr')

    fig.tight_layout(pad=1.0)
    fig.savefig('shell_grads.png',dpi=200)

    plt.show()
    


main()
