

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



def plot_render_rho(ax,time_Myr,sdf,xmin,xmax,ymin,ymax,vmin_in,vmax_in):
    '''
    Plot input sdf log col-dens on the given ax
    '''
    ax = sdf.render('rho', ax=ax, xlim=(xmin,xmax), ylim=(ymin,ymax), log_scale=True, cmap='gist_heat', vmin=vmin_in, vmax=vmax_in, cbar_kws=dict(label='column density [g/$\mathrm{cm}^2$]',orientation='vertical',shrink=1,pad=0.02))

    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.9*(ymax-ymin),r'$\mathrm{t_{SN}}$ = '+str(abs(round(time_Myr,3)))+' Myr', color='white')

    return 



def plot_render_u(ax,time_Myr,sdf,xmin,xmax,ymin,ymax,vmin_in,vmax_in):
    '''
    Plot input sdf log col-u on the given ax 
    '''
    ax = sdf.render('u', ax=ax, xlim=(xmin,xmax), ylim=(ymin,ymax), log_scale=True, cmap='gnuplot2', vmin=vmin_in, vmax=vmax_in, cbar_kws=dict(label=r'column internal energy [$\mathrm{erg/g\ cm}$]',orientation='vertical',shrink=1,pad=0.02), dens_weight=True)
    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.9*(ymax-ymin), r'$\mathrm{t_{SN}}$ = '+str(abs(round(time_Myr,3)))+' Myr', color='white')

    return 


fig = plt.figure(figsize=(10, 15),dpi=100)

rho_min = 1e-1
rho_max = 1e+4
u_min = 1e+10
u_max = 1e+14

xmin = -6
xmax = 6
ymin = -6
ymax = 6


time_Myr, sdf = read_dump('cloud_20_10_clrsink14_01350')
ax1 = fig.add_subplot(421)
plot_render_rho(ax1,time_Myr-t_SN,sdf,xmin,xmax,ymin,ymax,rho_min,rho_max)
ax2 = fig.add_subplot(422)
plot_render_u(ax2,time_Myr-t_SN,sdf,xmin,xmax,ymin,ymax,u_min,u_max)

time_Myr, sdf = read_dump('cloud_20_10_clrsink14_04450')
ax3 = fig.add_subplot(423)
plot_render_rho(ax3,time_Myr-t_SN,sdf,xmin,xmax,ymin,ymax,rho_min,rho_max)
ax4 = fig.add_subplot(424)
plot_render_u(ax4,time_Myr-t_SN,sdf,xmin,xmax,ymin,ymax,u_min,u_max)

time_Myr, sdf = read_dump('cloud_20_10_clrsink14_05000')
ax5 = fig.add_subplot(425)
plot_render_rho(ax5,time_Myr-t_SN,sdf,xmin,xmax,ymin,ymax,rho_min,rho_max)
ax6 = fig.add_subplot(426)
plot_render_u(ax6,time_Myr-t_SN,sdf,xmin,xmax,ymin,ymax,u_min,u_max)

time_Myr, sdf = read_dump('cloud_20_10_clrsink14_09570')
ax7 = fig.add_subplot(427)
plot_render_rho(ax7,time_Myr-t_SN,sdf,xmin,xmax,ymin,ymax,rho_min,rho_max)
ax8 = fig.add_subplot(428)
plot_render_u(ax8,time_Myr-t_SN,sdf,xmin,xmax,ymin,ymax,u_min,u_max)


fig.tight_layout(pad=0.1)

fig.savefig('cloud_density_u.png',dpi=200)
fig.savefig('cloud_density_u.pdf',format='pdf')

plt.show()

















