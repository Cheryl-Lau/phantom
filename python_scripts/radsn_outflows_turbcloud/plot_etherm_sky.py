
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

t_SN = 3.2573e-02*utime/Myr



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

    im = ax.pcolormesh(theta_mesh, phi_mesh, columnden, norm=colors.LogNorm(vmin=vmin_in,vmax=vmax_in), cmap='jet')
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, location='right')

    cbar.ax.set_ylabel(r'thermal energy $\mathrm{[erg\ sr^{-1}]}$', labelpad=15, rotation=270)
    ax.text(10,1.5,r'$\mathrm{t_{SN}}$ = '+str(abs(round(time_Myr,3)))+' Myr', fontsize=14)

    return 


fig = plt.figure(figsize=(15,15),dpi=200)

vmin = 3e13
vmax = 5e15


### 2 pc maps ###

ax1 = fig.add_subplot(421, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_01350')
plot_mollweide(ax1,'2pc_maps/sky_etherm_cloud_20_10_clrsink14_01350.dat',vmin,vmax,time_Myr-t_SN)
ax1.set_title('Within 2 pc radius', fontsize=18, pad=30)

ax2 = fig.add_subplot(423, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_04450')
plot_mollweide(ax2,'2pc_maps/sky_etherm_cloud_20_10_clrsink14_04450.dat',vmin,vmax,time_Myr-t_SN)

ax3 = fig.add_subplot(425, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_05000')
plot_mollweide(ax3,'2pc_maps/sky_etherm_cloud_20_10_clrsink14_05000.dat',vmin,vmax,time_Myr-t_SN)

ax4 = fig.add_subplot(427, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_09570')
plot_mollweide(ax4,'2pc_maps/sky_etherm_cloud_20_10_clrsink14_09570.dat',vmin,vmax,time_Myr-t_SN)


### 15 pc maps ###

ax5 = fig.add_subplot(422, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_01350')
plot_mollweide(ax5,'15pc_maps/sky_etherm_cloud_20_10_clrsink14_01350.dat',vmin,vmax,time_Myr-t_SN)
ax5.set_title('Within 15 pc radius', fontsize=18, pad=30)

ax6 = fig.add_subplot(424, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_04450')
plot_mollweide(ax6,'15pc_maps/sky_etherm_cloud_20_10_clrsink14_04450.dat',vmin,vmax,time_Myr-t_SN)

ax7 = fig.add_subplot(426, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_05000')
plot_mollweide(ax7,'15pc_maps/sky_etherm_cloud_20_10_clrsink14_05000.dat',vmin,vmax,time_Myr-t_SN)

ax8 = fig.add_subplot(428, projection='mollweide')
time_Myr, sdf = read_dump('cloud_20_10_clrsink14_09570')
plot_mollweide(ax8,'15pc_maps/sky_etherm_cloud_20_10_clrsink14_09570.dat',vmin,vmax,time_Myr-t_SN)



fig.tight_layout(pad=0.1)


plt.savefig('column_etherm_sky.png',dpi=200)
plt.savefig('column_etherm_sky.pdf',format='pdf')
plt.show()




























