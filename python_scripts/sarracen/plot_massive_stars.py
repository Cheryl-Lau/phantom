

import numpy as np
import matplotlib.pyplot as plt 
import sarracen


utime = 4.706e14
udist = 3.086e+19
umass = 1.989e+33
unit_density = 6.768e-23
unit_velocity = 6.558e+3
unit_ergg = 4.301e+7
Myr = 1e6*365*24*60*60 
solarm = 1.9891e+33


sdf, sdf_sinks = sarracen.read_phantom('cloud_21_20_11660')


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


def plot_render_rho_xy(ax,time_Myr,sdf,centre,radius):
    '''
    Plot input sdf log col-dens on the given ax
    '''
    ax = sdf.render('rho', ax=ax, xlim=(centre[0]-radius,centre[0]+radius), ylim=(centre[1]-radius,centre[1]+radius), log_scale=True, cmap='gist_heat', vmin=1e-1, vmax=1e5, cbar_kws=dict(label='column density [g$\ \mathrm{cm}^{-2}$]',orientation='vertical',shrink=0.7,pad=0.02))

    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(centre[0]-radius*0.8,centre[1]+radius*0.8,'t = '+str(round(time_Myr,3))+' Myr',color='white')

    return 


def plot_render_rho_xz(ax,time_Myr,sdf,centre,radius):
    '''
    Plot input sdf log col-dens on the given ax
    '''
    ax = sdf.render('rho', ax=ax, xlim=(centre[0]-radius,centre[0]+radius), ylim=(centre[2]-radius,centre[2]+radius), rotation=[0,0,90], rot_origin=[0,0,0], log_scale=True, cmap='gist_heat', vmin=1e-1, vmax=1e5, cbar_kws=dict(label='column density [g$\ \mathrm{cm}^{-2}$]',orientation='vertical',shrink=0.7,pad=0.02))

    ax.set_xlabel('x [pc]')
    ax.set_ylabel('z [pc]')
    ax.text(centre[0]-radius*0.8,centre[2]+radius*0.8,'t = '+str(round(time_Myr,3))+' Myr',color='white')

    return 


def plot_sinks_xy(ax,sdf_sinks,centre,radius): 

    ax.scatter(sdf_sinks['x'],sdf_sinks['y'],s=2,color='white')
    ax.set_xlim([centre[0]-radius,centre[0]+radius])
    ax.set_ylim((centre[1]-radius,centre[1]+radius))

    return 


def plot_sinks_xz(ax,sdf_sinks,centre,radius): 

    ax.scatter(sdf_sinks['x'],sdf_sinks['z'],s=2,color='white')
    ax.set_xlim([centre[0]-radius,centre[0]+radius])
    ax.set_ylim((centre[2]-radius,centre[2]+radius))

    return 


def plot_massive_sinks_xy(ax,sdf_sinks,centre,radius):

    m_threshold = 8.0
    imassive = np.where(sdf_sinks['m']*umass > m_threshold*solarm)[0]
    if (len(imassive) == 0):
        print('No massive stars found!!!')
        print('Biggest m: ',np.max(sdf_sinks['m']))
        return 

    x = sdf_sinks['x'].loc[imassive]
    y = sdf_sinks['y'].loc[imassive]
    m = sdf_sinks['m'].loc[imassive]

    m_max = np.max(m)*1.1

    pcm = ax.scatter(x,y,s=4,c=m,cmap='GnBu',marker='o',vmin=m_threshold,vmax=m_max)
    ax.set_xlim([centre[0]-radius,centre[0]+radius])
    ax.set_ylim((centre[1]-radius,centre[1]+radius))

    return pcm


def plot_massive_sinks_xz(ax,sdf_sinks,centre,radius):

    m_threshold = 8.0
    imassive = np.where(sdf_sinks['m']*umass > m_threshold*solarm)[0]
    if (len(imassive) == 0):
        print('No massive stars found!!!')
        print('Biggest m: ',np.max(sdf_sinks['m']))
        return 

    x = sdf_sinks['x'].loc[imassive]
    z = sdf_sinks['z'].loc[imassive]
    m = sdf_sinks['m'].loc[imassive]

    m_max = np.max(m)*1.1
    print(m_max)

    pcm = ax.scatter(x,z,s=4,c=m,cmap='GnBu',marker='o',vmin=m_threshold,vmax=m_max)
    ax.set_xlim([centre[0]-radius,centre[0]+radius])
    ax.set_ylim((centre[2]-radius,centre[2]+radius))

    return pcm



fig, axes = plt.subplots(ncols=2, subplot_kw=dict(frameon=True), figsize=(10,5))
plt.gca().set_aspect('equal')
ax1 = axes[0]
ax2 = axes[1]


t_Myr, sdf = process_sdf(sdf)
centre = [5,4,0]
radius = 15.0 

plot_render_rho_xy(ax1,t_Myr,sdf,centre,radius)
plot_render_rho_xz(ax2,t_Myr,sdf,centre,radius)

plot_sinks_xy(ax1,sdf_sinks,centre,radius)
plot_sinks_xz(ax2,sdf_sinks,centre,radius)
pcm1 = plot_massive_sinks_xy(ax1,sdf_sinks,centre,radius)
pcm2 = plot_massive_sinks_xz(ax2,sdf_sinks,centre,radius)

cbar1 = fig.colorbar(pcm1, ax=ax1, location='top', shrink=0.75)
cbar1.set_label(r'sink mass [$M_\odot$]')

cbar2 = fig.colorbar(pcm2, ax=ax2, location='top', shrink=0.75)
cbar2.set_label(r'sink mass [$M_\odot$]')


fig.tight_layout(pad=1.0)


plt.savefig('massive_stars.png',dpi=200)
plt.show()

















