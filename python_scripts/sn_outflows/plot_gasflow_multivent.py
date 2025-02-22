'''
Plotting the velocities at vent for varying number of channels 
'''

import numpy as np 
import matplotlib.pyplot as plt 
import glob
import os
import math

utime = 4.706e14
udist = 3.086e+19
umass = 1.989e+33
unit_density = 6.768e-23
unit_velocity = 6.558e+3
unit_ergg = 4.301e+7
Myr = 1e6*365*24*60*60


# Shock arrival time
t0 = 0.092*Myr 



def sort_array(array_ref,array_in):
    return [i for _, i in sorted(zip(array_ref, array_in))]


def plot_analytical_model(ax,labels,colours):

    time_modch1,vx_modch1,ptherm_modch1,pram_modch1 = np.loadtxt('semi_confined_model.txt',unpack=True)
    time_modch2,vx_modch2,ptherm_modch2,pram_modch2 = np.loadtxt('semi_confined_2chnl_model.txt',unpack=True)
    time_modch3,vx_modch3,ptherm_modch3,pram_modch3 = np.loadtxt('semi_confined_3chnl_model.txt',unpack=True)
    time_modch4,vx_modch4,ptherm_modch4,pram_modch4 = np.loadtxt('semi_confined_4chnl_model.txt',unpack=True)

    time_modch1 = time_modch1/Myr
    time_modch2 = time_modch2/Myr
    time_modch3 = time_modch3/Myr
    time_modch4 = time_modch4/Myr 

    ax.plot(time_modch1,vx_modch1,color='black',linewidth=1,linestyle='dotted',label=labels[0]+' model')
    ax.plot(time_modch2,vx_modch2,color='black',linewidth=1,linestyle='dashed',label=labels[1]+' model')
    ax.plot(time_modch3,vx_modch3,color='black',linewidth=1,linestyle='dashdot',label=labels[2]+' model')
    ax.plot(time_modch4,vx_modch4,color='black',linewidth=1,linestyle='solid',label=labels[3]+' model')


    time_modff,vx_modff,ptherm_modff,pram_modff = np.loadtxt('free_field_model.txt',unpack=True)
    time_modff = time_modff/Myr 
    
    ax.plot(time_modff,vx_modff,color='blue',linewidth=1,linestyle='dashed',label='free-field model')

    return 


def plot_sim_output(ax,chnl_dirs,labels,colours): 

    c = 0 
    for chnl_dir in chnl_dirs:

        time_all = []
        vx_all = []
        ptherm_all = []
        pram_all = []

        path = chnl_dir+"/gasflow_semiconf_*.dat"
        for filename in glob.glob(path):  

            time_cgs,x,y,z = np.loadtxt(filename,max_rows=1,skiprows=1,unpack=True)
            rho,vx,u,ptherm,pram = np.loadtxt(filename,max_rows=3,skiprows=3,unpack=True)

            time_all.append((time_cgs-t0)/Myr)
            vx_all.append(vx)
            ptherm_all.append(ptherm)
            pram_all.append(pram)

        time_sorted = sort_array(time_all,time_all)
        vx_sorted = sort_array(time_all,vx_all)
        ptherm_sorted = sort_array(time_all,ptherm_all)
        pram_sorted = sort_array(time_all,pram_all)

        ax.scatter(time_sorted,vx_sorted,s=2,color=colours[c],label=labels[c])
        c += 1 

    return 



fig,(ax1,ax2) = plt.subplots(nrows=2,sharex=False,subplot_kw=dict(frameon=True),figsize=(5,7),dpi=200)


## Vary number of channels ##
chnl_dirs = ['semi_confined/vary_chnlsize/chnl005','semi_confined/vary_nchnl/2chnl','semi_confined/vary_nchnl/3chnl','semi_confined/vary_nchnl/4chnl']

labels = ['1 channel','2 channels','3 channels','4 channels']
colours = ['springgreen','orange','lightskyblue','lightcoral']

plot_analytical_model(ax1,labels,colours)
plot_sim_output(ax1,chnl_dirs,labels,colours)

ax1.set_xlim([-0.01,0.38])
ax1.set_ylim([4E5,9E8])
ax1.set_yscale('log')
ax1.set_xlabel('time [Myr]')
ax1.set_ylabel('velocity [cm $\mathrm{s^{-1}}$]')
ax1.legend(fontsize=7.5,loc='upper right')


## Vary channel size ##
chnl_dirs = ['semi_confined/vary_chnlsize/chnl005','semi_confined/vary_chnlsize/chnl010','semi_confined/vary_chnlsize/chnl015','semi_confined/vary_chnlsize/chnl020']

labels = [r'$\Omega = 0.05 \times 4 \pi\ \mathrm{sr}$',r'$\Omega = 0.10 \times 4 \pi\ \mathrm{sr}$',r'$\Omega = 0.15 \times 4 \pi\ \mathrm{sr}$',r'$\Omega = 0.20 \times 4 \pi\ \mathrm{sr}$']
colours = ['indigo','firebrick','orange','yellowgreen']

plot_analytical_model(ax2,labels,colours)
plot_sim_output(ax2,chnl_dirs,labels,colours)

ax2.set_xlim([-0.01,0.38])
ax2.set_ylim([4E5,9E8])
ax2.set_yscale('log')
ax2.set_xlabel('time [Myr]')
ax2.set_ylabel('velocity [cm $\mathrm{s^{-1}}$]')
ax2.legend(fontsize=7.5,loc='upper right')



fig.tight_layout(pad=1.0)

fig.savefig('vel_multivent_adiabatic.pdf',format='pdf')
fig.savefig('vel_multivent_adiabatic.png')





