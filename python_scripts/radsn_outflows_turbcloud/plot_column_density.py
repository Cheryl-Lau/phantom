
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import sarracen


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

    # Convert units 
    time = time*utime
    time_Myr = time/(1e6*365*24*60*60)
    sdf['rho'] = sdf['rho']*unit_density 
    sdf['u'] = sdf['u']*unit_ergg

    return time_Myr, sdf



def plot_mollweide(ax,thetaphi_file,vmin_in,vmax_in,time_Myr):

    columnden = np.loadtxt(thetaphi_file)

    theta = np.linspace(-np.pi, np.pi,np.shape(columnden)[1])
    phi = np.linspace(-np.pi/2., np.pi/2.,np.shape(columnden)[0])
    theta_mesh, phi_mesh = np.meshgrid(theta,phi)

    im = ax.pcolormesh(theta_mesh, phi_mesh, columnden, norm=colors.LogNorm(vmin=vmin_in,vmax=vmax_in), cmap='inferno')
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, location='right')

    cbar.ax.set_ylabel(r'column density $\mathrm{[g\ cm^{-2}]}$', labelpad=15, rotation=270)
    ax.set_title('t = '+str(round(time_Myr,3))+' Myr',loc='right')

    return 


fig = plt.figure(figsize=(15,12),dpi=200)

vmin = 5e-5
vmax = 1e-2 

ax1 = fig.add_subplot(321, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_01350')
plot_mollweide(ax1,'sky_density_cloud_20_10_clrsink14_01350.dat',vmin,vmax,time_Myr)

ax2 = fig.add_subplot(323, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_04450')
plot_mollweide(ax2,'sky_density_cloud_20_10_clrsink14_04450.dat',vmin,vmax,time_Myr)

ax3 = fig.add_subplot(325, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_05000')
plot_mollweide(ax3,'sky_density_cloud_20_10_clrsink14_05000.dat',vmin,vmax,time_Myr)

ax4 = fig.add_subplot(322, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_06000')
plot_mollweide(ax4,'sky_density_cloud_20_10_clrsink14_06000.dat',vmin,vmax,time_Myr)

ax5 = fig.add_subplot(324, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_08000')
plot_mollweide(ax5,'sky_density_cloud_20_10_clrsink14_08000.dat',vmin,vmax,time_Myr)

ax6 = fig.add_subplot(326, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_09570')
plot_mollweide(ax6,'sky_density_cloud_20_10_clrsink14_09570.dat',vmin,vmax,time_Myr)



fig.tight_layout(pad=0.1)


plt.savefig('column_density.png',dpi=200)
plt.savefig('column_density.pdf',format='pdf')
plt.show()






















