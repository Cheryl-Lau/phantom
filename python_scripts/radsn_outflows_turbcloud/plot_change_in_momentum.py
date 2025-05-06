
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
solarm = 1.989e+33


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
    sdf['vx'] = sdf['vx']*unit_velocity
    sdf['vy'] = sdf['vy']*unit_velocity
    sdf['vz'] = sdf['vz']*unit_velocity
    sdf['vmag'] = sdf['vmag']*unit_velocity

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
    loc = [x,y,z]
  
    return loc


def extract_highrho(sdf): 

    sdf = sdf.loc[sdf['rho'] > 1e-24]

    return sdf




def plot_change_in_mu(ax_r,ax_rho,sn_loc,t_sn,sdf1,sdf2):

    sdf1_new = sdf1.loc[sdf1.index.isin(sdf2.index)]
    sdf2 = sdf2.loc[sdf2.index.isin(sdf1.index)]
    sdf1 = sdf1_new 

    x1 = sdf1['x'] - sn_loc[0]
    y1 = sdf1['y'] - sn_loc[1]
    z1 = sdf1['z'] - sn_loc[2]
    sdf1['r'] = np.sqrt(x1**2 + y1**2 + z1**2) 
    x_unitvec = x1/sdf1['r']
    y_unitvec = y1/sdf1['r']
    z_unitvec = z1/sdf1['r']
    sdf1['vr'] = sdf1['vx']*x_unitvec + sdf1['vy']*y_unitvec + sdf1['vz']*z_unitvec

    x2 = sdf2['x'] - sn_loc[0]
    y2 = sdf2['y'] - sn_loc[1]
    z2 = sdf2['z'] - sn_loc[2]
    sdf2['r'] = np.sqrt(x2**2 + y2**2 + z2**2) 
    x_unitvec = x2/sdf2['r']
    y_unitvec = y2/sdf2['r']
    z_unitvec = z2/sdf2['r']
    sdf2['vr'] = sdf2['vx']*x_unitvec + sdf2['vy']*y_unitvec + sdf2['vz']*z_unitvec

    pmass = 1e-2*solarm
    dmu = pmass*(sdf2['vr']-sdf1['vr'])


    # Plot mu against r 

    pcm = ax_r.scatter(sdf1['r'], dmu, c=sdf2['T'],norm=matplotlib.colors.LogNorm(vmin=1e1,vmax=1e3), s=1, alpha=0.5)

    ax_r.set_xlim([1.5e-2,4e+1])
    ax_r.set_ylim([1e+30,1e+42])

    ax_r.set_xlabel('$r_0$ [pc]')
    ax_r.set_ylabel(r'$m(v_r-v_{r,0})\ \mathrm{[g\ cm\ s^{-1}]}$')
    ax_r.set_xscale('log')
    ax_r.set_yscale('log')
    ax_r.text(4e-2, 1e+41, r'$t_\mathrm{SN} = $'+str(round(t_sn,2))+' Myr', color='black')


    # Plot mu against rho

    pcm = ax_rho.scatter(sdf1['rho'], dmu, c=sdf2['T'],norm=matplotlib.colors.LogNorm(vmin=1e1,vmax=1e3), s=1, alpha=0.5)

    ax_rho.set_xlim([1e-26,1e-14])
    ax_rho.set_ylim([1e+30,1e+42])

    ax_rho.set_xlabel(r'$\rho_0\ \mathrm{[g\ cm^{-3}]}$')
    ax_rho.set_ylabel(r'$m(v_r-v_{r,0})\ \mathrm{[g\ cm\ s^{-1}]}$')
    ax_rho.set_xscale('log')
    ax_rho.set_yscale('log')
    ax_rho.text(1e-25, 1e+41, r'$t_\mathrm{SN} = $'+str(round(t_sn,2))+' Myr', color='black')


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



fig, axes = plt.subplots(nrows=2, ncols=3, subplot_kw=dict(frameon=True), figsize=(10,5.5))


# Semi-confined SN 

ax_r = axes[0,0]
ax_rho = axes[1,0]

time1, sdf1, sdf1_sinks = read_dump('turbcloud_semiconf/cloud_20_10_clrsink14_01350')
sdf1 = limit_box(sdf1)
time2, sdf2, sdf2_sinks = read_dump('turbcloud_semiconf/cloud_20_10_clrsink14_08000')
t_sn = time2 - time1
sn_loc = get_sink_loc(14,sdf1_sinks)

pcm = plot_change_in_mu(ax_r,ax_rho,sn_loc,t_sn,sdf1,sdf2)


cbar = fig.colorbar(pcm, ax=ax_r, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)
cbar = fig.colorbar(pcm, ax=ax_rho, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)

ax_r.set_title('Semi-confined SN (incl. envelope)')



# Turb FF - smear 

ax_r = axes[0,1]
ax_rho = axes[1,1]

time1, sdf1, sdf1_sinks = read_dump('turbff_smear/turbffsmear_00000')
sdf1 = limit_box(sdf1)
time2, sdf2, sdf2_sinks = read_dump('turbff_smear/turbffsmear_08000')
t_sn = time2 - time1
sn_loc = [0,0,0]

pcm = plot_change_in_mu(ax_r,ax_rho,sn_loc,t_sn,sdf1,sdf2)


cbar = fig.colorbar(pcm, ax=ax_r, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)
cbar = fig.colorbar(pcm, ax=ax_rho, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)

ax_r.set_title(r'Turb FF - $\rho_\mathrm{smear}$')



# Turb FF - env 

ax_r = axes[0,2]
ax_rho = axes[1,2]

time1, sdf1, sdf1_sinks = read_dump('turbff_env2/turbffenv_00000')
sdf1 = limit_box(sdf1)
time2, sdf2, sdf2_sinks = read_dump('turbff_env2/turbffenv_00200')
t_sn = time2 - time1
sn_loc = [0,0,0]

pcm = plot_change_in_mu(ax_r,ax_rho,sn_loc,t_sn,sdf1,sdf2)


cbar = fig.colorbar(pcm, ax=ax_r, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)
cbar = fig.colorbar(pcm, ax=ax_rho, cmap='viridis', location='right', shrink=0.9,pad=0.02)
cbar.set_label('Temperature [K]',rotation=270,labelpad=15)

ax_r.set_title(r'Turb FF - $\rho_\mathrm{env}$')




fig.tight_layout()
fig.subplots_adjust(wspace=0.45)
fig.subplots_adjust(hspace=0.45)

plt.savefig('momentum_change.png',dpi=200)
plt.show()














