
import numpy as np
import matplotlib.pyplot as plt 
import sarracen
import inspect
import glob

import matplotlib
matplotlib.use("AGG")


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
    ax = sdf.render('rho', ax=ax, xlim=(xmin,xmax), ylim=(ymin,ymax), log_scale=True, cmap='gist_heat', vmin=vmin_in, vmax=vmax_in, cbar_kws=dict(label='column density [g/$\mathrm{cm}^2$]',orientation='vertical',shrink=0.6,pad=0.02))

    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.9*(ymax-ymin),r'$\mathrm{t_{SN}}$ = '+str('{:.3f}'.format(round(abs(time_Myr),3)))+' Myr', color='white')

    return 



def plot_render_u(ax,time_Myr,sdf,xmin,xmax,ymin,ymax,vmin_in,vmax_in):
    '''
    Plot input sdf log col-u on the given ax 
    '''
    ax = sdf.render('u', ax=ax, xlim=(xmin,xmax), ylim=(ymin,ymax), log_scale=True, cmap='gnuplot2', vmin=vmin_in, vmax=vmax_in, cbar_kws=dict(label=r'column internal energy [$\mathrm{erg/g\ cm}$]',orientation='vertical',shrink=0.6,pad=0.02), dens_weight=True)
    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.9*(ymax-ymin), r'$\mathrm{t_{SN}}$ = '+str('{:.3f}'.format(round(abs(time_Myr),3)))+' Myr', color='white')

    return 



rho_min = 1e-1
rho_max = 1e+4
u_min = 1e+10
u_max = 1e+14

xmin = -7
xmax = 7
ymin = -6
ymax = 6


for dumpfile in sorted(glob.glob('./cloud_20_10_clrsink14_*')):  
    print(dumpfile)

    fig = plt.figure(figsize=(10, 4),dpi=200)

    time_Myr, sdf = read_dump(dumpfile)
    ax1 = fig.add_subplot(121)
    plot_render_rho(ax1,time_Myr-t_SN,sdf,xmin,xmax,ymin,ymax,rho_min,rho_max)
    ax2 = fig.add_subplot(122)
    plot_render_u(ax2,time_Myr-t_SN,sdf,xmin,xmax,ymin,ymax,u_min,u_max)

    fig.tight_layout(pad=0.1)

    imfilename = 'render_'+dumpfile[-5:]+'.png'
    print(imfilename)
    fig.savefig(imfilename,dpi=200)















