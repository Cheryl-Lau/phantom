
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
    iselect = np.where(np.logical_and(df['radius'] == rad,df['chnl']==0.1))[0]
    time_cf01 = df['time'].iloc[iselect].to_numpy()
    etot_cf01 = df['etot'].iloc[iselect].to_numpy()
    #iselect = np.where(np.logical_and(df['radius'] == rad,df['chnl']==0.2))[0]
    #time_cf02 = df['time'].iloc[iselect].to_numpy()
    #etot_cf02 = df['etot'].iloc[iselect].to_numpy()
    iselect = np.where(np.logical_and(df['radius'] == rad,df['chnl']==0.3))[0]
    time_cf03 = df['time'].iloc[iselect].to_numpy()
    etot_cf03 = df['etot'].iloc[iselect].to_numpy()

    time0_ff = np.min(time_ff00)   # sync to time when shock arrives rad 
    time0_cf = np.min(time_cf01)  

    ax.scatter(time_ff00-time0_ff,etot_ff00,s=1,color='darkblue',label='free-field')
    ax.scatter(time_cf01-time0_cf,etot_cf01,s=1,color='salmon',label='confined '+'$\Omega = 0.1 \cdot 4 \pi$')
    #ax.scatter(time_cf02-time0_cf,etot_cf02,s=1,color='orangered',label='confined '+'$\Omega = 0.2 \cdot 4 \pi$')
    ax.scatter(time_cf03-time0_cf,etot_cf03,s=1,color='firebrick',label='confined '+'$\Omega = 0.3 \cdot 4 \pi$')
    ax.set_yscale('log')
    ax.set_xlabel('time [Myr]')
    ax.set_xlim([0,3.1])
    ax.set_ylim([3e46,6e51])
    ax.legend(loc='lower left',fontsize=9)
    ax.text(1.8,1.5e51,'r < '+str(int(rad))+' pc')

    ax.minorticks_on()
    ax.xaxis.set_tick_params(which='minor')

    return 


df = read_input_combineframe()

fig, (ax1,ax2,ax3,ax4) = plt.subplots(ncols=4, sharey=True, subplot_kw=dict(frameon=True), figsize=(10,3), dpi=300)

plot_energ_in_rad(ax1,20.0,df)
plot_energ_in_rad(ax2,50.0,df)
plot_energ_in_rad(ax3,80.0,df)
plot_energ_in_rad(ax4,100.0,df)

ax1.set_ylabel('total energy [erg]')
fig.tight_layout(pad=1.0)


plt.savefig('energy_insphere.png')











