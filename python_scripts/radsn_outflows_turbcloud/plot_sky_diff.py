
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



def plot_mollweide_diff(ax,thetaphi_file1,thetaphi_file2,vmin_in,vmax_in):

    columnden1 = np.loadtxt(thetaphi_file1)
    columnden1 = columnden1/umass

    columnden2 = np.loadtxt(thetaphi_file2)
    columnden2 = columnden2/umass

    theta = np.linspace(-np.pi, np.pi,np.shape(columnden1)[1])
    phi = np.linspace(-np.pi/2., np.pi/2.,np.shape(columnden1)[0])
    theta_mesh, phi_mesh = np.meshgrid(theta,phi)

    # log difference 
#    columnden_logdiff = np.log10(columnden1/columnden2)
#    columnden_diff = 10**columnden_logdiff
    columnden_diff = (columnden1 - columnden2)/columnden1

    im = ax.pcolormesh(theta_mesh, phi_mesh, columnden_diff, norm=colors.LogNorm(vmin=vmin_in,vmax=vmax_in), cmap='cividis')
    cbar = fig.colorbar(im, ax=ax, shrink=0.5, location='right')

    cbar.ax.set_ylabel(r'fractional mass removed $\mathrm{[M_\odot\ sr^{-1}]}$', labelpad=15, rotation=270)

    return 


fig = plt.figure(figsize=(15,5),dpi=200)

vmin = 2e-1  # 1.6e0
vmax = 1e0   #  2.7e0 


### 2 pc maps ###

ax1 = fig.add_subplot(121, projection='mollweide')
map1 = '2pc_maps/sky_density_cloud_20_10_clrsink14_01350.dat'
map2 = '2pc_maps/sky_density_cloud_20_10_clrsink14_04450.dat'
plot_mollweide_diff(ax1,map1,map2,vmin,vmax)
ax1.set_title('Within 2 pc radius', fontsize=18, pad=30)


### 15 pc maps ###

ax2 = fig.add_subplot(122, projection='mollweide')
map1 = '15pc_maps/sky_density_cloud_20_10_clrsink14_01350.dat'
map2 = '15pc_maps/sky_density_cloud_20_10_clrsink14_04450.dat'
plot_mollweide_diff(ax2,map1,map2,vmin,vmax)
ax2.set_title('Within 15 pc radius', fontsize=18, pad=30)




fig.tight_layout(pad=0.1)


plt.savefig('diff_sky.png',dpi=200)
plt.savefig('diff_sky.pdf',format='pdf')
plt.show()























