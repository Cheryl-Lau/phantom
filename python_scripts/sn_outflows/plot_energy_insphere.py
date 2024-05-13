
'''
Plots the evolution of total energy contained within spheres of different radii
'''

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import glob 
import os 

df = pd.DataFrame(columns=('chnl', 'time', 'radius', 'ekin','etherm','etot'))

Myr = 1E6*365*24*60*60 
pc = 3.086E+18

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

#for filename in glob.glob('semi_confined02/energy_insphere_withhii_chnl02_*.dat'): 
#    time = np.loadtxt(filename,skiprows=1,max_rows=1,unpack=True)
#    radius, ke, te, tote = np.loadtxt(filename,skiprows=3,unpack=True)
#    df02 = pd.DataFrame({'chnl':[0.2]*len(radius),'time':time/Myr,'radius':radius/pc,'ekin':ke,'etherm':te,'etot':tote})
#    df = pd.concat([df,df02],ignore_index=True)

for filename in glob.glob('semi_confined03/energy_insphere_withhii_chnl03_*.dat'): 
    time = np.loadtxt(filename,skiprows=1,max_rows=1,unpack=True)
    radius, ke, te, tote = np.loadtxt(filename,skiprows=3,unpack=True)
    df03 = pd.DataFrame({'chnl':[0.3]*len(radius),'time':time/Myr,'radius':radius/pc,'ekin':ke,'etherm':te,'etot':tote})
    df = pd.concat([df,df03],ignore_index=True)

print(df)



fig, (ax1,ax2,ax3) = plt.subplots(nrows=3, sharex=True, subplot_kw=dict(frameon=True), figsize=(6,10), dpi=300)
plt.subplots_adjust(hspace=.0)


# rad = 20 pc
iselect = np.where(np.logical_and(df['radius'] == 20.0,df['chnl']==0.0))[0]
time_r20_ff00 = df['time'].iloc[iselect].to_numpy()
etot_r20_ff00 = df['etot'].iloc[iselect].to_numpy()
iselect = np.where(np.logical_and(df['radius'] == 20.0,df['chnl']==0.1))[0]
time_r20_cf01 = df['time'].iloc[iselect].to_numpy()
etot_r20_cf01 = df['etot'].iloc[iselect].to_numpy()
#iselect = np.where(np.logical_and(df['radius'] == 20.0,df['chnl']==0.2))[0]
#time_r20_cf02 = df['time'].iloc[iselect].to_numpy()
#etot_r20_cf02 = df['etot'].iloc[iselect].to_numpy()
iselect = np.where(np.logical_and(df['radius'] == 20.0,df['chnl']==0.3))[0]
time_r20_cf03 = df['time'].iloc[iselect].to_numpy()
etot_r20_cf03 = df['etot'].iloc[iselect].to_numpy()

time0_r20_ff = np.min(time_r20_ff00)   # sync to time when shock arrives rad 
time0_r20_cf = np.min(time_r20_cf01)  

ax1 = fig.add_subplot(311)
ax1.scatter(time_r20_ff00-time0_r20_ff,etot_r20_ff00,s=1,color='darkblue',label='free-field')
ax1.scatter(time_r20_cf01-time0_r20_cf,etot_r20_cf01,s=1,color='salmon',label='$\Omega = 0.1 \cdot 4 \pi$')
#ax1.scatter(time_r20_cf02-time0_r20_cf,etot_r20_cf02,s=1,color='orangered',label='$\Omega = 0.2 \cdot 4 \pi$')
ax1.scatter(time_r20_cf03-time0_r20_cf,etot_r20_cf03,s=1,color='firebrick',label='$\Omega = 0.3 \cdot 4 \pi$')
ax1.set_yscale('log')
ax1.set_ylim([3E49,1E52])
ax1.set_xlabel('Time [Myr]')
ax1.set_ylabel('Total energy [erg]')
ax1.text(0.1,5E51,'radius = 20 pc')
ax1.legend()


# rad = 50 pc
iselect = np.where(np.logical_and(df['radius'] == 50.0,df['chnl']==0.0))[0]
time_r50_ff00 = df['time'].iloc[iselect].to_numpy()
etot_r50_ff00 = df['etot'].iloc[iselect].to_numpy()
iselect = np.where(np.logical_and(df['radius'] == 50.0,df['chnl']==0.1))[0]
time_r50_cf01 = df['time'].iloc[iselect].to_numpy()
etot_r50_cf01 = df['etot'].iloc[iselect].to_numpy()
#iselect = np.where(np.logical_and(df['radius'] == 50.0,df['chnl']==0.2))[0]
#time_r50_cf02 = df['time'].iloc[iselect].to_numpy()
#etot_r50_cf02 = df['etot'].iloc[iselect].to_numpy()
iselect = np.where(np.logical_and(df['radius'] == 50.0,df['chnl']==0.3))[0]
time_r50_cf03 = df['time'].iloc[iselect].to_numpy()
etot_r50_cf03 = df['etot'].iloc[iselect].to_numpy()

time0_r50_ff = np.min(time_r50_ff00)   # sync to time when shock arrives rad 
time0_r50_cf = np.min(time_r50_cf01)  

ax2 = fig.add_subplot(312)
ax2.scatter(time_r50_ff00-time0_r50_ff,etot_r50_ff00,s=1,color='darkblue',label='free-field')
ax2.scatter(time_r50_cf01-time0_r50_cf,etot_r50_cf01,s=1,color='salmon',label='$\Omega = 0.1 \cdot 4 \pi$')
#ax2.scatter(time_r50_cf02-time0_r50_cf,etot_r50_cf02,s=1,color='orangered',label='$\Omega = 0.2 \cdot 4 \pi$')
ax2.scatter(time_r50_cf03-time0_r50_cf,etot_r50_cf03,s=1,color='firebrick',label='$\Omega = 0.3 \cdot 4 \pi$')
ax2.set_yscale('log')
ax2.set_ylim([1E50,5E52])
ax2.set_xlabel('Time [Myr]')
ax2.set_ylabel('Total energy [erg]')
ax2.text(0.1,2.5E52,'radius = 50 pc')
ax2.legend()


# rad = 100 pc
iselect = np.where(np.logical_and(df['radius'] == 100.0,df['chnl']==0.0))[0]
time_r100_ff00 = df['time'].iloc[iselect].to_numpy()
etot_r100_ff00 = df['etot'].iloc[iselect].to_numpy()
iselect = np.where(np.logical_and(df['radius'] == 100.0,df['chnl']==0.1))[0]
time_r100_cf01 = df['time'].iloc[iselect].to_numpy()
etot_r100_cf01 = df['etot'].iloc[iselect].to_numpy()
#iselect = np.where(np.logical_and(df['radius'] == 100.0,df['chnl']==0.2))[0]
#time_r100_cf02 = df['time'].iloc[iselect].to_numpy()
#etot_r100_cf02 = df['etot'].iloc[iselect].to_numpy()
iselect = np.where(np.logical_and(df['radius'] == 100.0,df['chnl']==0.3))[0]
time_r100_cf03 = df['time'].iloc[iselect].to_numpy()
etot_r100_cf03 = df['etot'].iloc[iselect].to_numpy()

time0_r100_ff = np.min(time_r100_ff00)   # sync to time when shock arrives rad 
time0_r100_cf = np.min(time_r100_cf01)  

ax3 = fig.add_subplot(313)
ax3.scatter(time_r100_ff00-time0_r100_ff,etot_r100_ff00,s=1,color='darkblue',label='free-field')
ax3.scatter(time_r100_cf01-time0_r100_cf,etot_r100_cf01,s=1,color='salmon',label='$\Omega = 0.1 \cdot 4 \pi$')
#ax3.scatter(time_r100_cf02-time0_r100_cf,etot_r100_cf02,s=1,color='orangered',label='$\Omega = 0.2 \cdot 4 \pi$')
ax3.scatter(time_r100_cf03-time0_r100_cf,etot_r100_cf03,s=1,color='firebrick',label='$\Omega = 0.3 \cdot 4 \pi$')
ax3.set_yscale('log')
ax3.set_ylim([3E50,5E52])
ax3.set_xlabel('Time [Myr]')
ax3.set_ylabel('Total energy [erg]')
ax3.text(0.1,2.5E52,'radius = 100 pc')
ax3.legend()


plt.savefig('energy_insphere.png')
plt.show()












