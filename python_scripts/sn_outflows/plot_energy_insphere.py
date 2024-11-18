
'''
Plots the evolution of total energy contained within spheres of different radii
'''

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import glob 
import os 


Myr = 1E6*365*24*60*60 
pc = 3.086E+18


def plot_line(ax,x,y,colour):
    '''
    For connecting the scatter dots 
    sorts x, then sort y according to x 
    '''
    x_sorted = sorted(x)
    y_sorted = [ydum for _, ydum in sorted(zip(x,y))]
    x_sorted.pop(0)
    y_sorted.pop(0)
    ax.plot(x_sorted,y_sorted,color=colour,linewidth=1)    

    return 


def plot_dashed_line(ax,x,y,colour):
    '''
    For total energy  
    sorts x, then sort y according to x 
    '''
    x_sorted = sorted(x)
    y_sorted = [ydum for _, ydum in sorted(zip(x,y))]
    x_sorted.pop(0)
    y_sorted.pop(0)

    ax.plot(x_sorted,y_sorted,'--',color=colour,linewidth=1)    

    return 


def read_input_combineframe():

    df = pd.DataFrame(columns=('chnl', 'time', 'radius', 'ekin','etherm','etot'))

    for filename in glob.glob('free_field/energy_insphere_nocloud_bigenv_*.dat'):  # over_time 
        time = np.loadtxt(filename,skiprows=1,max_rows=1,unpack=True)
        radius, ke, te, tote = np.loadtxt(filename,skiprows=3,unpack=True)
        df00 = pd.DataFrame({'chnl':[0.0]*len(radius),'time':time/Myr,'radius':radius/pc,'ekin':ke,'etherm':te,'etot':tote})
        df = pd.concat([df,df00],ignore_index=True)

    for filename in glob.glob('semi_confined01/energy_insphere_withhii_chnl01_*.dat'): 
        time = np.loadtxt(filename,skiprows=1,max_rows=1,unpack=True)
        radius, ke, te, tote = np.loadtxt(filename,skiprows=3,unpack=True)
        df01 = pd.DataFrame({'chnl':[0.1]*len(radius),'time':time/Myr,'radius':radius/pc,'ekin':ke,'etherm':te,'etot':tote})
        df = pd.concat([df,df01],ignore_index=True)

    for filename in glob.glob('semi_confined02/energy_insphere_withhii_chnl02_*.dat'): 
        time = np.loadtxt(filename,skiprows=1,max_rows=1,unpack=True)
        radius, ke, te, tote = np.loadtxt(filename,skiprows=3,unpack=True)
        df02 = pd.DataFrame({'chnl':[0.2]*len(radius),'time':time/Myr,'radius':radius/pc,'ekin':ke,'etherm':te,'etot':tote})
        df = pd.concat([df,df02],ignore_index=True)

    for filename in glob.glob('semi_confined03/energy_insphere_withhii_chnl03_*.dat'): 
        time = np.loadtxt(filename,skiprows=1,max_rows=1,unpack=True)
        radius, ke, te, tote = np.loadtxt(filename,skiprows=3,unpack=True)
        df03 = pd.DataFrame({'chnl':[0.3]*len(radius),'time':time/Myr,'radius':radius/pc,'ekin':ke,'etherm':te,'etot':tote})
        df = pd.concat([df,df03],ignore_index=True)

    print(df)

    return df


def plot_energ_in_rad(ax,rad,df):

    iselect = np.where(np.logical_and(df['radius'] == rad,df['chnl']==0.0))[0]
    time_ff00 = df['time'].iloc[iselect].to_numpy()
    etot_ff00 = df['etot'].iloc[iselect].to_numpy()
    ekin_ff00 = df['ekin'].iloc[iselect].to_numpy()
    etherm_ff00 = df['etherm'].iloc[iselect].to_numpy()
    iselect = np.where(np.logical_and(df['radius'] == rad,df['chnl']==0.1))[0]
    time_cf01 = df['time'].iloc[iselect].to_numpy()
    etot_cf01 = df['etot'].iloc[iselect].to_numpy()
    ekin_cf01 = df['ekin'].iloc[iselect].to_numpy()
    etherm_cf01 = df['etherm'].iloc[iselect].to_numpy()
    iselect = np.where(np.logical_and(df['radius'] == rad,df['chnl']==0.3))[0]
    time_cf03 = df['time'].iloc[iselect].to_numpy()
    etot_cf03 = df['etot'].iloc[iselect].to_numpy()
    ekin_cf03 = df['ekin'].iloc[iselect].to_numpy()
    etherm_cf03 = df['etherm'].iloc[iselect].to_numpy()

    time0_ff = np.min(time_ff00)   # sync to time when shock arrives rad 
    time0_cf = np.min(time_cf01)  


    ax.scatter(time_ff00-time0_ff,ekin_ff00,s=1,color='navy',label='free-field $\mathrm{E_{kin}}$')
    ax.scatter(time_cf01-time0_cf,ekin_cf01,s=1,color='orangered',label='confined $\mathrm{E_{kin}}$ '+'$\Omega = 0.1 \cdot 4 \pi$')
    ax.scatter(time_cf03-time0_cf,ekin_cf03,s=1,color='orange',label='confined $\mathrm{E_{kin}}$ '+'$\Omega = 0.3 \cdot 4 \pi$')

    ax.scatter(time_ff00-time0_ff,etherm_ff00,s=1,color='cornflowerblue',label='free-field $\mathrm{E_{therm}}$')
    ax.scatter(time_cf01-time0_cf,etherm_cf01,s=1,color='coral',label='confined $\mathrm{E_{therm}}$ '+'$\Omega = 0.1 \cdot 4 \pi$')
    ax.scatter(time_cf03-time0_cf,etherm_cf03,s=1,color='gold',label='confined $\mathrm{E_{therm}}$ '+'$\Omega = 0.3 \cdot 4 \pi$')

    plot_line(ax,time_ff00-time0_ff,ekin_ff00,'navy')
    plot_line(ax,time_ff00-time0_ff,etherm_ff00,'cornflowerblue')
    plot_line(ax,time_cf01-time0_cf,ekin_cf01,'orangered')
    plot_line(ax,time_cf01-time0_cf,etherm_cf01,'coral')
    plot_line(ax,time_cf03-time0_cf,ekin_cf03,'orange')
    plot_line(ax,time_cf03-time0_cf,etherm_cf03,'gold')

    # total energy 
    plot_dashed_line(ax,time_ff00-time0_ff,etot_ff00,'midnightblue')
    plot_dashed_line(ax,time_cf01-time0_cf,etot_cf01,'darkred')

    ax.set_yscale('log')
    ax.set_xlabel('time [Myr]')
    ax.set_xlim([0,3.1])
    ax.set_ylim([8e45,9e51])
    ax.text(1.4,2e51,'radius < '+str(int(rad))+' pc')

    ax.minorticks_on()
    ax.xaxis.set_tick_params(which='minor')

    return 


df = read_input_combineframe()

fig, (ax1,ax2,ax3,ax4) = plt.subplots(ncols=4, sharey=True, subplot_kw=dict(frameon=True), figsize=(10,3), dpi=200)
#plt.subplots_adjust(hspace=.0)

plot_energ_in_rad(ax1,20.0,df)
plot_energ_in_rad(ax2,50.0,df)
plot_energ_in_rad(ax3,80.0,df)
plot_energ_in_rad(ax4,100.0,df)

ax1.set_ylabel('energy [erg]')
ax4.legend(loc='lower left',fontsize=7)

fig.tight_layout(pad=0.5)


plt.savefig('energy_insphere.png')











