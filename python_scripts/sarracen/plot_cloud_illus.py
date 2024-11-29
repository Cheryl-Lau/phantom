

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


sdf, sdf_sinks = sarracen.read_phantom('cloud_20_10_clsink139_v2_00600')


def process_sdf(sdf):
    '''
    Calculate extra quantities
    Convert to useful units 
    '''
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


def plot_render_rho(ax,time_Myr,sdf,centre,radius):
    '''
    Plot input sdf log col-dens on the given ax
    '''
    ax = sdf.render('rho', ax=ax, xlim=(centre[0]-radius,centre[0]+radius), ylim=(centre[1]-radius,centre[1]+radius), log_scale=True, cmap='gist_heat', vmin=1e1, vmax=1e5, cbar_kws=dict(label='column density [g/$\mathrm{cm}^2$]',orientation='vertical',shrink=0.7,pad=0.02))

    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(centre[0]-radius*0.8,centre[1]+radius*0.8,'t = '+str(round(time_Myr,3))+' Myr',color='white')

    return 


def plot_render_u(ax,time_Myr,sdf,centre,radius):
    '''
    Plot input sdf log col-u on the given ax 
    '''
    ax = sdf.render('u', ax=ax, xlim=(centre[0]-radius,centre[0]+radius), ylim=(centre[1]-radius,centre[1]+radius), log_scale=True, cmap='gnuplot2', vmin=1e10, vmax=1e14, cbar_kws=dict(label='column internal energy [erg/g cm]',orientation='vertical',shrink=0.7,pad=0.02), dens_weight=True)

    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(centre[0]-radius*0.8,centre[1]+radius*0.8,'t = '+str(round(time_Myr,3))+' Myr',color='white')

    return 


def get_sink_info(isink,sdf_sinks):

    irow = isink - 1 
    x = sdf_sinks['x'].iloc[irow]
    y = sdf_sinks['y'].iloc[irow]
    z = sdf_sinks['z'].iloc[irow]
    loc_sink = [x,y,z]
    mass_sink = sdf_sinks['m'].iloc[irow]
    hacc_sink = sdf_sinks['h'].iloc[irow]
    
    return loc_sink, mass_sink, hacc_sink


def plot_sinks(ax,sdf_sinks,centre,radius): 

    ax.scatter(sdf_sinks['x'], sdf_sinks['y'], s=2, color='white')
    ax.set_xlim([centre[0]-radius,centre[0]+radius])
    ax.set_ylim((centre[1]-radius,centre[1]+radius))

    return 


fig, (ax1,ax2) = plt.subplots(ncols=2, sharey=False, subplot_kw=dict(frameon=True), figsize=(10,5))

plt.gca().set_aspect('equal')


t_Myr, sdf = process_sdf(sdf)
centre, mass, hacc = get_sink_info(139,sdf_sinks)

print(centre,mass)
radius = 4.0 

plot_render_rho(ax1,t_Myr,sdf,centre,radius)
plot_render_u(ax2,t_Myr,sdf,centre,radius)
plot_sinks(ax1,sdf_sinks,centre,radius)
plot_sinks(ax2,sdf_sinks,centre,radius)

fig.tight_layout(pad=1.0)


plt.savefig('cloud_hiiregion_illus.png',dpi=200)
plt.show()
















