

import numpy as np
import matplotlib.pyplot as plt 
import sarracen
import inspect


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



def plot_render_rho(ax,time_Myr,sdf,xmin,xmax,ymin,ymax,vmin_in,vmax_in):
    '''
    Plot input sdf log col-dens on the given ax
    '''
    ax = sdf.render('rho', ax=ax, xlim=(xmin,xmax), ylim=(ymin,ymax), log_scale=True, cmap='gist_heat', vmin=vmin_in, vmax=vmax_in, cbar_kws=dict(label='column density [g/$\mathrm{cm}^2$]',orientation='vertical',shrink=0.85,pad=0.02))

    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.9*(ymax-ymin),'t = '+str(round(time_Myr,3))+' Myr', color='white')

    return 



def plot_render_u(ax,time_Myr,sdf,xmin,xmax,ymin,ymax,vmin_in,vmax_in):
    '''
    Plot input sdf log col-u on the given ax 
    '''
    ax = sdf.render('u', ax=ax, xlim=(xmin,xmax), ylim=(ymin,ymax), log_scale=True, cmap='gnuplot2', vmin=vmin_in, vmax=vmax_in, cbar_kws=dict(label=r'column internal energy [$\mathrm{erg/g\ cm}$]',orientation='vertical',shrink=0.85,pad=0.02), dens_weight=True)
    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.9*(ymax-ymin), 't = '+str(round(time_Myr,3))+' Myr', color='white')

    return 


fig1 = plt.figure(figsize=(9, 12),dpi=100)
fig2 = plt.figure(figsize=(9, 12),dpi=100)

rho_min = 1e-1
rho_max = 1e+4
u_min = 1e+10
u_max = 1e+14

xmin = -6
xmax = 6
ymin = -6
ymax = 6


time_Myr, sdf = read_dump('cloud_20_10_clrsink14_01350')
ax11 = fig1.add_subplot(321)
plot_render_rho(ax11,time_Myr,sdf,xmin,xmax,ymin,ymax,rho_min,rho_max)
ax21 = fig2.add_subplot(321)
plot_render_u(ax21,time_Myr,sdf,xmin,xmax,ymin,ymax,u_min,u_max)

time_Myr, sdf = read_dump('cloud_20_10_clrsink14_04450')
ax13 = fig1.add_subplot(323)
plot_render_rho(ax13,time_Myr,sdf,xmin,xmax,ymin,ymax,rho_min,rho_max)
ax23 = fig2.add_subplot(323)
plot_render_u(ax23,time_Myr,sdf,xmin,xmax,ymin,ymax,u_min,u_max)

time_Myr, sdf = read_dump('cloud_20_10_clrsink14_05000')
ax15 = fig1.add_subplot(325)
plot_render_rho(ax15,time_Myr,sdf,xmin,xmax,ymin,ymax,rho_min,rho_max)
ax25 = fig2.add_subplot(325)
plot_render_u(ax25,time_Myr,sdf,xmin,xmax,ymin,ymax,u_min,u_max)

time_Myr, sdf = read_dump('cloud_20_10_clrsink14_06000')
ax12 = fig1.add_subplot(322)
plot_render_rho(ax12,time_Myr,sdf,xmin,xmax,ymin,ymax,rho_min,rho_max)
ax22 = fig2.add_subplot(322)
plot_render_u(ax22,time_Myr,sdf,xmin,xmax,ymin,ymax,u_min,u_max)

time_Myr, sdf = read_dump('cloud_20_10_clrsink14_08000')
ax14 = fig1.add_subplot(324)
plot_render_rho(ax14,time_Myr,sdf,xmin,xmax,ymin,ymax,rho_min,rho_max)
ax24 = fig2.add_subplot(324)
plot_render_u(ax24,time_Myr,sdf,xmin,xmax,ymin,ymax,u_min,u_max)

time_Myr, sdf = read_dump('cloud_20_10_clrsink14_09570')
ax16 = fig1.add_subplot(326)
plot_render_rho(ax16,time_Myr,sdf,xmin,xmax,ymin,ymax,rho_min,rho_max)
ax26 = fig2.add_subplot(326)
plot_render_u(ax26,time_Myr,sdf,xmin,xmax,ymin,ymax,u_min,u_max)


fig1.tight_layout(pad=0.1)
fig2.tight_layout(pad=0.1)


fig1.savefig('cloud_density.png',dpi=200)
fig1.savefig('cloud_density.pdf',format='pdf')
fig2.savefig('cloud_internalenergy.png',dpi=200)
fig2.savefig('cloud_internalenergy.pdf',format='pdf')
plt.show()

















