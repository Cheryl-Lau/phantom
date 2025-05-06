
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sarracen
import pandas as pd
import matplotlib.gridspec as gridspec
from matplotlib._layoutgrid import plot_children


utime = 4.706e14
udist = 3.086e+19
umass = 1.989e+33
unit_density = 6.768e-23
unit_velocity = 6.558e+3
unit_ergg = 4.301e+7
Myr = 1e6*365*24*60*60 



def read_dump(dumpfile):
    '''
    Calculate extra quantities
    Convert to useful units 
    '''
    sdf, sdf_sinks = sarracen.read_phantom(dumpfile)

    # Get dumpfile info 
    params_dict = sdf.params 
    time = params_dict['time']
    npart = params_dict['npartoftype']
    gamma = params_dict['gamma']
    pmass = params_dict['massoftype']

    # Calculate extra quantities 
    sdf.calc_density()
    sdf['vmag'] = np.sqrt(sdf['vx']**2 + sdf['vy']**2 + sdf['vz']**2)
    sdf['P'] = sdf['u'] * sdf['rho'] * (sdf.params['gamma'] - 1.0)
    sdf['T'] = get_temp(sdf['u'])

    # Convert units 
    time = time*utime
    time_Myr = time/(1e6*365*24*60*60)
    sdf['rho'] = sdf['rho']*unit_density 
    sdf['u'] = sdf['u']*unit_ergg

    return time_Myr, sdf, sdf_sinks



def get_temp(u):

    kboltz = 1.38066e-16
    gmw = 2.381
    gamma = 5.0/3.0
    mass_proton_cgs = 1.67262158e-24

    return u/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg



def get_sink_loc(isink,sdf_sinks):

    x = sdf_sinks['x'].loc[isink]
    y = sdf_sinks['y'].loc[isink]
    z = sdf_sinks['z'].loc[isink]
    m = sdf_sinks['m'].loc[isink]
    sn_loc = [x,y,z]
  
    return sn_loc


def extract_highrho(sdf): 

    sdf = sdf.loc[sdf['rho'] > 1e-24]

    return sdf

def extract_hightemp(sdf):

    sdf = sdf.loc[sdf['T'] > 1e4]

    return sdf




def plot_change_in_r(ax,sn_loc,t_sn,sdf1,sdf2):

    sdf1_new = sdf1.loc[sdf1.index.isin(sdf2.index)]
    sdf2 = sdf2.loc[sdf2.index.isin(sdf1.index)]
    sdf1 = sdf1_new 

    x1 = sdf1['x'] - sn_loc[0]
    y1 = sdf1['y'] - sn_loc[1]
    z1 = sdf1['z'] - sn_loc[2]
    sdf1['r'] = np.sqrt(x1**2 + y1**2 + z1**2) 

    x2 = sdf2['x'] - sn_loc[0]
    y2 = sdf2['y'] - sn_loc[1]
    z2 = sdf2['z'] - sn_loc[2]
    sdf2['r'] = np.sqrt(x2**2 + y2**2 + z2**2) 

    sdf1['dr'] = abs(sdf2['r'] - sdf1['r'])

    pcm = ax.scatter(sdf1['r'], sdf1['dr'], c=sdf2['T'],norm=matplotlib.colors.LogNorm(vmin=1e1,vmax=1e3), s=1, alpha=0.5)

    ax.set_xlim([2e-2,2.5e+1])
    ax.set_ylim([1e-8,1e+3])

    ax.set_xlabel('$r_0$ [pc]')
    ax.set_ylabel('$|r-r_0|$ [pc]')
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.text(3e-2, 5e+1, r'$t_\mathrm{SN} = $'+str(round(t_sn,2))+' Myr', color='black')

    return pcm




def plot_change_in_rho(ax,sn_loc,t_sn,sdf1,sdf2):

    sdf1_new = sdf1.loc[sdf1.index.isin(sdf2.index)]
    sdf2 = sdf2.loc[sdf2.index.isin(sdf1.index)]
    sdf1 = sdf1_new 

    x1 = sdf1['x'] - sn_loc[0]
    y1 = sdf1['y'] - sn_loc[1]
    z1 = sdf1['z'] - sn_loc[2]
    sdf1['r'] = np.sqrt(x1**2 + y1**2 + z1**2) 

    sdf1['drho'] = sdf2['rho'] #/sdf1['rho']

    pcm = ax.scatter(sdf1['rho'], sdf1['drho'], c=sdf2['T'], norm=matplotlib.colors.LogNorm(vmin=1e1,vmax=1e3), s=1, alpha=0.5)

    ax.set_xlim([1e-26,1e-15])
    ax.set_ylim([1e-30,1e-15])

    rholine = np.logspace(-26,-15,100)
    ax.plot(rholine,rholine,color='black',linewidth=1,label=r'$\rho = \rho_0$')

    ax.set_xlabel(r'$\rho_0\ \mathrm{[g\ cm^{-3}]}$')
    ax.set_ylabel(r'$\rho\ \mathrm{[g\ cm^{-3}]}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='lower right')

    ax.text(5e-26, 1e-17, r'$t_\mathrm{SN} = $'+str(round(t_sn,2))+' Myr', color='black')

    return pcm


def limit_box(sdf):

    xmin = -10
    xmax = 10 
    ymin = -10
    ymax = 10
    zmin = -10 
    zmax = 10
    iinbox = np.where((sdf['x']>xmin) & (sdf['x']<xmax) & (sdf['y']>ymin) & (sdf['y']<ymax) & (sdf['z']>zmin) & (sdf['z']<zmax))[0]

    sdf = sdf.iloc[iinbox]

    return sdf




fig = plt.figure(figsize=(10,10),layout="constrained")
gs = fig.add_gridspec(4,3)


# Semi-confined SN 

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])

time1, sdf1, sdf1_sinks = read_dump('turbcloud_semiconf/cloud_20_10_clrsink14_01350')
sdf1 = limit_box(sdf1)
time2, sdf2, sdf2_sinks = read_dump('turbcloud_semiconf/cloud_20_10_clrsink14_04500')
t_sn = time2 - time1
sn_loc = get_sink_loc(14,sdf1_sinks)

pcm = plot_change_in_r(ax1,sn_loc,t_sn,sdf1,sdf2)
pcm = plot_change_in_rho(ax2,sn_loc,t_sn,sdf1,sdf2)

cbar = fig.colorbar(pcm, ax=ax1, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)
cbar = fig.colorbar(pcm, ax=ax2, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)

ax1.set_title('Semi-confined SN (incl. envelope)')



# FF - cloud

ax3 = fig.add_subplot(gs[0, 1])
ax4 = fig.add_subplot(gs[1, 1])

time1, sdf1, sdf1_sinks = read_dump('ff_cloud/ffcloud_00000')
sdf1 = limit_box(sdf1)
time2, sdf2, sdf2_sinks = read_dump('ff_cloud/ffcloud_02000')
t_sn = time2 - time1
sn_loc = [0,0,0]

pcm = plot_change_in_r(ax3,sn_loc,t_sn,sdf1,sdf2)
pcm = plot_change_in_rho(ax4,sn_loc,t_sn,sdf1,sdf2)

cbar = fig.colorbar(pcm, ax=ax3, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)
cbar = fig.colorbar(pcm, ax=ax4, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)

ax3.set_title(r'FF - $\rho_\mathrm{cloud}$')



# Turb FF - cloud 

ax5 = fig.add_subplot(gs[0, 2])
ax6 = fig.add_subplot(gs[1, 2])

time1, sdf1, sdf1_sinks = read_dump('turbff_cloud/turbffcloud_00000')
sdf1 = limit_box(sdf1)
time2, sdf2, sdf2_sinks = read_dump('turbff_cloud/turbffcloud_01900')
t_sn = time2 - time1
sn_loc = [0,0,0]

pcm = plot_change_in_r(ax5,sn_loc,t_sn,sdf1,sdf2)
pcm = plot_change_in_rho(ax6,sn_loc,t_sn,sdf1,sdf2)

cbar = fig.colorbar(pcm, ax=ax5, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)
cbar = fig.colorbar(pcm, ax=ax6, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)

ax5.set_title(r'Turb FF - $\rho_\mathrm{cloud}$')



# FF - env

ax7 = fig.add_subplot(gs[2, 1])
ax8 = fig.add_subplot(gs[3, 1])

time1, sdf1, sdf1_sinks = read_dump('ff_env3/freefield_00000')
sdf1 = limit_box(sdf1)
time2, sdf2, sdf2_sinks = read_dump('ff_env3/freefield_00200')
t_sn = time2 - time1
sn_loc = [0,0,0]

pcm = plot_change_in_r(ax7,sn_loc,t_sn,sdf1,sdf2)
pcm = plot_change_in_rho(ax8,sn_loc,t_sn,sdf1,sdf2)

cbar = fig.colorbar(pcm, ax=ax7, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)
cbar = fig.colorbar(pcm, ax=ax8, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)

ax7.set_title(r'FF - $\rho_\mathrm{env}$')



# Turb FF - env 

ax9 = fig.add_subplot(gs[2, 2])
ax10 = fig.add_subplot(gs[3, 2])

time1, sdf1, sdf1_sinks = read_dump('turbff_env2/turbffenv_00000')
sdf1 = limit_box(sdf1)
time2, sdf2, sdf2_sinks = read_dump('turbff_env2/turbffenv_00200')
t_sn = time2 - time1
sn_loc = [0,0,0]

pcm = plot_change_in_r(ax9,sn_loc,t_sn,sdf1,sdf2)
pcm = plot_change_in_rho(ax10,sn_loc,t_sn,sdf1,sdf2)

cbar = fig.colorbar(pcm, ax=ax9, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)
cbar = fig.colorbar(pcm, ax=ax10, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)

ax9.set_title(r'Turb FF - $\rho_\mathrm{env}$')



plt.savefig('sn_kicks_vs_cloud_env.png',dpi=200)















