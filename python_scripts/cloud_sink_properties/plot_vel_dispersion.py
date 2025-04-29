

import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
from scipy.optimize import curve_fit
import sarracen


utime = 4.706e14
udist = 3.086e+19
umass = 1.989e+33
unit_density = 6.768e-23
unit_velocity = 6.558e+3
unit_ergg = 4.301e+7
Myr = 1e6*365*24*60*60 
solarm = 1.9891e+33



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

    # Convert units 
    time = time*utime
    time_Myr = time/(1e6*365*24*60*60)
    sdf['rho'] = sdf['rho']*unit_density 
    sdf['u'] = sdf['u']*unit_ergg

    sdf_sinks = sdf_sinks[sdf_sinks['m'] > 0]

    return time_Myr, sdf, sdf_sinks



def plot_vel_dispersion(ax,veldisfile):

    radii = np.loadtxt(veldisfile,max_rows=1)
    data = np.loadtxt(veldisfile,skiprows=1,usecols=np.arange(1,len(radii)+1))

    df = pd.DataFrame(data, columns=radii.tolist())
    df = df.iloc[:, np.arange(10,len(radii))]
    df = df[df!=0]

    radii = df.columns.tolist()
    sigma = df.mean(axis=0)

    ax.plot(radii,sigma,color='black')

    fit_larson(ax,radii,sigma)

    return 


def plot_settings(ax,rho0str,time_Myr,t_ff):

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([3.5e-2,1.5e+1])
    ax.set_ylim([9e3,4e5])
    ax.set_xlabel('length-scale $R$')
    ax.set_ylabel('velocity dispersion $\sigma$')
    
    ax.text(4.5e-2,2.8e+5,rho0str+r'$\ \ \ \alpha =$'+str(format(alpha,'.2f')),color='black')
    ax.text(4.5e-2,2.0e+5,'t = '+str(round(time_Myr,3))+' Myr ('+str(round(time_Myr/t_ff,3))+' $t_{ff}$)',color='black')

    ax.grid(linestyle='--',color='darkgrey',linewidth=0.4,which='minor')
    ax.grid(linestyle='--',color='darkgrey',linewidth=0.7,which='major')

    return 


def fit_larson(ax,radii,sigmas):

    r = np.logspace(np.log10(min(radii)),np.log10(max(radii)),num=100)
    maxpt = np.where(np.array(radii) < 1e0)[0][-1]
    param, covar = curve_fit(larson_relation, radii[0:maxpt], sigmas.tolist()[0:maxpt])
    C = param[0]

    ax.plot(r,larson_relation(r, C),'--',color='grey',label=r'Larson relation $\sigma \propto R^{1/2}$')

    return 


def larson_relation(r, C):

    sigma = C * r**(0.5)

    return sigma 



fig, axes = plt.subplots(nrows=3, ncols=3, subplot_kw=dict(frameon=True), figsize=(10,10))


##### rho 1e-22 #####

rho0str = r'$\rho_0 = 10^{-22}\ \mathrm{g\ cm^{-3}}$'
t_ff = 6.662 # Myr 

ax22_07 = axes[0,0]
alpha = 0.7
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_07/cloud_22_07_07025')
plot_vel_dispersion(ax22_07,'elongev_22_07/velocity_dispersion_cloud_22_07_07025.dat')
plot_settings(ax22_07,rho0str,time_Myr,t_ff)

ax22_10 = axes[1,0]
alpha = 1.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_10/cloud_22_10_08100')
plot_vel_dispersion(ax22_10,'elongev_22_10/velocity_dispersion_cloud_22_10_08000.dat')
plot_settings(ax22_10,rho0str,time_Myr,t_ff)

ax22_20 = axes[2,0]
alpha = 2.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_20/cloud_22_20_08700')
plot_vel_dispersion(ax22_20,'elongev_22_20/velocity_dispersion_cloud_22_20_08700.dat')
plot_settings(ax22_20,rho0str,time_Myr,t_ff)


##### rho 1e-21 #####

rho0str = r'$\rho_0 = 10^{-21}\ \mathrm{g\ cm^{-3}}$'
t_ff = 2.107 # Myr 

ax21_07 = axes[0,1]
alpha = 0.7
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_07/cloud_21_07_02462')
plot_vel_dispersion(ax21_07,'elongev_21_07/velocity_dispersion_cloud_21_07_02462.dat')
plot_settings(ax21_07,rho0str,time_Myr,t_ff)

ax21_10 = axes[1,1]
alpha = 1.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_10/cloud_21_10_12010')
plot_vel_dispersion(ax21_10,'elongev_21_10/velocity_dispersion_cloud_21_10_12010.dat')
plot_settings(ax21_10,rho0str,time_Myr,t_ff)

ax21_20 = axes[2,1]
alpha = 2.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_20/cloud_21_20_10900')
plot_vel_dispersion(ax21_20,'elongev_21_20/velocity_dispersion_cloud_21_20_10900.dat')
plot_settings(ax21_20,rho0str,time_Myr,t_ff)


##### rho 1e-20 #####

rho0str = r'$\rho_0 = 10^{-20}\ \mathrm{g\ cm^{-3}}$'
t_ff = 0.666 # Myr 

ax20_07 = axes[0,2]
alpha = 0.7
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_07/cloud_20_07_00585')
plot_vel_dispersion(ax20_07,'elongev_20_07/velocity_dispersion_cloud_20_07_00585.dat')
plot_settings(ax20_07,rho0str,time_Myr,t_ff)
ax20_07.legend(loc='lower right')

ax20_10 = axes[1,2]
alpha = 1.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_10/cloud_20_10_03748')
plot_vel_dispersion(ax20_10,'elongev_20_10/velocity_dispersion_cloud_20_10_03748.dat')
plot_settings(ax20_10,rho0str,time_Myr,t_ff)
ax20_10.legend(loc='lower right')

ax20_20 = axes[2,2]
alpha = 2.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_20/cloud_20_20_02800')
plot_vel_dispersion(ax20_20,'elongev_20_20/velocity_dispersion_cloud_20_20_02800.dat')
plot_settings(ax20_20,rho0str,time_Myr,t_ff)
ax20_20.legend(loc='lower right')


fig.tight_layout()


plt.savefig('clouds_vel_dispersion.png',dpi=200)
plt.savefig('clouds_vel_dispersion.pdf',format='pdf')

plt.show()
















