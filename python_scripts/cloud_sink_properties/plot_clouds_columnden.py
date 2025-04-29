
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

    return time_Myr, sdf, sdf_sinks




def plot_render_rho(ax,time_Myr,sdf,t_ff,alpha,rho0str,xmin,xmax,ymin,ymax,rho_min,rho_max):
    '''
    Plot input sdf log col-dens on the given ax
    '''
    ax = sdf.render('rho', ax=ax, xlim=(xmin,xmax), ylim=(ymin,ymax), log_scale=True, cmap='gist_heat', vmin=rho_min, vmax=rho_max, cbar_kws=dict(label='',orientation='vertical',shrink=0.5,pad=0.02))

    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]') 
    ax.text(xmin+(xmax-xmin)*0.1,ymin+(ymax-ymin)*0.9,'t = '+str(round(time_Myr,3))+' Myr ('+str(round(time_Myr/t_ff,3))+' $t_{ff}$)',color='white')
    ax.text(xmin+(xmax-xmin)*0.05,ymin+(ymax-ymin)*0.05,rho0str+r'$\ \ \ \alpha =$'+str(format(alpha,'.2f')),color='white')

    return 



def plot_sinks(ax,sdf_sinks,xmin,xmax,ymin,ymax): 

    ax.scatter(sdf_sinks['x'],sdf_sinks['y'],s=2,color='white')
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])

    return 



def plot_massive_sinks(ax,sdf_sinks,xmin,xmax,ymin,ymax):

    m_threshold = 8.0
    imassive = np.where(sdf_sinks['m']*umass > m_threshold*solarm)[0]
    if (len(imassive) == 0):
        print('No massive stars found!!!')
        print('Biggest m: ',np.max(sdf_sinks['m']))
        return 

    x = sdf_sinks['x'].loc[imassive]
    y = sdf_sinks['y'].loc[imassive]
    m = sdf_sinks['m'].loc[imassive]

    m_max = np.max(m)
    imax = sdf_sinks.index[sdf_sinks['m'] == m_max] 
    
    isink_max = imax[0] + 1 
    x_max = sdf_sinks['x'].loc[imax]
    y_max = sdf_sinks['y'].loc[imax]

    print('isink:',isink_max)
    print('x',x_max.to_numpy())
    print('y',y_max.to_numpy())
    print('z',sdf_sinks['z'].loc[imax].to_numpy())
    print('m',sdf_sinks['m'].loc[imax].to_numpy())

    pcm = ax.scatter(x,y,s=4,c=m,cmap='winter',marker='o',vmin=m_threshold,vmax=m_max)
    ax.scatter(x_max,y_max,s=40,c=m_max,cmap='winter',marker='*',vmin=m_threshold,vmax=m_max)

    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])

    return pcm



fig, axes = plt.subplots(nrows=3, ncols=3, subplot_kw=dict(frameon=True), figsize=(12,12))
plt.gca().set_aspect('equal')

rho_min = 1e-1
rho_max = 1e5


##### rho 1e-22 #####

rho0str = r'$\rho_0 = 10^{-22}\ \mathrm{g\ cm^{-3}}$'
t_ff = 6.662 # Myr 

ax22_07 = axes[0,0]
alpha = 0.7 
xmin = -2
xmax = 12 
ymin = 2
ymax = 16
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_07/cloud_22_07_07025')
plot_render_rho(ax22_07,time_Myr,sdf,t_ff,alpha,rho0str,xmin,xmax,ymin,ymax,rho_min,rho_max)
plot_sinks(ax22_07,sdf_sinks,xmin,xmax,ymin,ymax)
pcm = plot_massive_sinks(ax22_07,sdf_sinks,xmin,xmax,ymin,ymax)
cbar = fig.colorbar(pcm, ax=ax22_07, location='top', shrink=0.6)
cbar.set_label(r'sink mass [$M_\odot$]')

ax22_10 = axes[1,0]
alpha = 1.0 
xmin = -2
xmax = 12 
ymin = 5
ymax = 19
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_10/cloud_22_10_08100')
plot_render_rho(ax22_10,time_Myr,sdf,t_ff,alpha,rho0str,xmin,xmax,ymin,ymax,rho_min,rho_max)
plot_sinks(ax22_10,sdf_sinks,xmin,xmax,ymin,ymax)
pcm = plot_massive_sinks(ax22_10,sdf_sinks,xmin,xmax,ymin,ymax)
cbar = fig.colorbar(pcm, ax=ax22_10, location='top', shrink=0.6)
cbar.set_label(r'sink mass [$M_\odot$]')

ax22_20 = axes[2,0]
alpha = 2.0 
xmin = -6
xmax = 20 
ymin = 9
ymax = 35
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_20/cloud_22_20_08700')
plot_render_rho(ax22_20,time_Myr,sdf,t_ff,alpha,rho0str,xmin,xmax,ymin,ymax,rho_min,rho_max)
plot_sinks(ax22_20,sdf_sinks,xmin,xmax,ymin,ymax)
pcm = plot_massive_sinks(ax22_20,sdf_sinks,xmin,xmax,ymin,ymax)
cbar = fig.colorbar(pcm, ax=ax22_20, location='top', shrink=0.6)
cbar.set_label(r'sink mass [$M_\odot$]')



##### rho 1e-21 #####

rho0str = r'$\rho_0 = 10^{-21}\ \mathrm{g\ cm^{-3}}$'
t_ff = 2.107 # Myr 

ax21_07 = axes[0,1]
alpha = 0.7 
xmin = -6
xmax = 6 
ymin = -4
ymax = 8
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_07/cloud_21_07_02462')
plot_render_rho(ax21_07,time_Myr,sdf,t_ff,alpha,rho0str,xmin,xmax,ymin,ymax,rho_min,rho_max)
plot_sinks(ax21_07,sdf_sinks,xmin,xmax,ymin,ymax)
pcm = plot_massive_sinks(ax21_07,sdf_sinks,xmin,xmax,ymin,ymax)
cbar = fig.colorbar(pcm, ax=ax21_07, location='top', shrink=0.6)
cbar.set_label(r'sink mass [$M_\odot$]')

ax21_10 = axes[1,1]
alpha = 1.0
xmin = -4
xmax = 8 
ymin = -3
ymax = 9
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_10/cloud_21_10_12010')
plot_render_rho(ax21_10,time_Myr,sdf,t_ff,alpha,rho0str,xmin,xmax,ymin,ymax,rho_min,rho_max)
plot_sinks(ax21_10,sdf_sinks,xmin,xmax,ymin,ymax)
pcm = plot_massive_sinks(ax21_10,sdf_sinks,xmin,xmax,ymin,ymax)
cbar = fig.colorbar(pcm, ax=ax21_10, location='top', shrink=0.6)
cbar.set_label(r'sink mass [$M_\odot$]')

ax21_20 = axes[2,1]
alpha = 2.0
xmin = -8
xmax = 16
ymin = -4
ymax = 20
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_20/cloud_21_20_10900')
plot_render_rho(ax21_20,time_Myr,sdf,t_ff,alpha,rho0str,xmin,xmax,ymin,ymax,rho_min,rho_max)
plot_sinks(ax21_20,sdf_sinks,xmin,xmax,ymin,ymax)
pcm = plot_massive_sinks(ax21_20,sdf_sinks,xmin,xmax,ymin,ymax)
cbar = fig.colorbar(pcm, ax=ax21_20, location='top', shrink=0.6)
cbar.set_label(r'sink mass [$M_\odot$]')



##### rho 1e-20 #####

rho0str = r'$\rho_0 = 10^{-20}\ \mathrm{g\ cm^{-3}}$'
t_ff = 0.666 # Myr 

ax20_07 = axes[0,2]
alpha = 0.7 
xmin = -4
xmax = 4
ymin = -3
ymax = 5
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_07/cloud_20_07_00585')
plot_render_rho(ax20_07,time_Myr,sdf,t_ff,alpha,rho0str,xmin,xmax,ymin,ymax,rho_min,rho_max)
plot_sinks(ax20_07,sdf_sinks,xmin,xmax,ymin,ymax)
pcm = plot_massive_sinks(ax20_07,sdf_sinks,xmin,xmax,ymin,ymax)
cbar = fig.colorbar(pcm, ax=ax20_07, location='top', shrink=0.6)
cbar.set_label(r'sink mass [$M_\odot$]')

ax20_10 = axes[1,2]
alpha = 1.0
xmin = -4
xmax = 4
ymin = -3
ymax = 5
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_10/cloud_20_10_03748')
plot_render_rho(ax20_10,time_Myr,sdf,t_ff,alpha,rho0str,xmin,xmax,ymin,ymax,rho_min,rho_max)
plot_sinks(ax20_10,sdf_sinks,xmin,xmax,ymin,ymax)
pcm = plot_massive_sinks(ax20_10,sdf_sinks,xmin,xmax,ymin,ymax)
cbar = fig.colorbar(pcm, ax=ax20_10, location='top', shrink=0.6)
cbar.set_label(r'sink mass [$M_\odot$]')

ax20_20 = axes[2,2]
alpha = 2.0
xmin = -4
xmax = 4
ymin = -3
ymax = 5
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_20/cloud_20_20_02800')
plot_render_rho(ax20_20,time_Myr,sdf,t_ff,alpha,rho0str,xmin,xmax,ymin,ymax,rho_min,rho_max)
plot_sinks(ax20_20,sdf_sinks,xmin,xmax,ymin,ymax)
pcm = plot_massive_sinks(ax20_20,sdf_sinks,xmin,xmax,ymin,ymax)
cbar = fig.colorbar(pcm, ax=ax20_20, location='top', shrink=0.6)
cbar.set_label(r'sink mass [$M_\odot$]')



fig.tight_layout()
fig.subplots_adjust(wspace=0.12)
fig.subplots_adjust(hspace=0.12)

plt.savefig('clouds_columnden.png',dpi=200)
plt.savefig('clouds_columnden.pdf',format='pdf')
plt.show()










