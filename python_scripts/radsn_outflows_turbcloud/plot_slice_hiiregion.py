
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



def plot_render_u(ax,time_Myr,sdf,xmin,xmax,ymin,ymax,zmin,zmax,vmin_in,vmax_in):
    '''
    Plot input sdf log col-u on the given ax 
    '''
    ax = sdf.render('u', ax=ax, xlim=(xmin,xmax), ylim=(ymin,ymax), log_scale=True, cmap='gnuplot2', vmin=vmin_in, vmax=vmax_in, cbar_kws=dict(label=r'column internal energy [$\mathrm{erg/g\ cm}$]',orientation='vertical',shrink=0.9,pad=0.02), dens_weight=False)
    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.9*(ymax-ymin), 't = '+str(round(time_Myr,3))+' Myr', color='white')
    ax.text(xmax - 0.4*(xmax-xmin), ymin + 0.9*(ymax-ymin), str(zmin)+' pc < z < '+str(zmax)+' pc', color='white')

    return 



def plot_slice_nH(ax,time_Myr,xmin,xmax,ymin,ymax,zmin,zmax):

    max_m = max(m)
    m_scaled = m/max_m

    islice = np.where(np.logical_and(z > zmin,z < zmax))
    x_slice = x[islice]
    y_slice = y[islice]
    m_slice = m[islice]
    h_slice = h[islice]
    nH_slice = nH[islice]

    max_m = max(m_slice)
    m_scaled = m_slice/max_m

    pcm = ax.scatter(x_slice,y_slice,s=1,c=nH_slice,alpha=0.6)
        
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.9*(ymax-ymin),'time = '+str(round(time_Myr,3))+' Myr')
    ax.text(xmax - 0.4*(xmax-xmin), ymin + 0.9*(ymax-ymin),str(zmin)+' pc < z < '+str(zmax)+' pc')

    # Mark the star 
    ax.scatter(centre[0],centre[1],marker='*',s=70,color='salmon',label='ionizing source')
    ax.legend(loc='lower right')

    return pcm



fig, axes = plt.subplots(2, 2, sharex=False, sharey=False, subplot_kw=dict(frameon=True), figsize=(12,10))
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]

radius = 4.5
centre = [7.349882733E-01, 9.668074170E-01, 9.430135498E-02] # sink 14
xmin = centre[0]-radius
xmax = centre[0]+radius
ymin = centre[1]-radius
ymax = centre[1]+radius

#u_min = 1e+10  # density-weighted
#u_max = 1e+14
u_min = 1e+30
u_max = 1e+33

zmin1 = -0.7
zmax1 = -0.3

zmin2 = -0.3
zmax2 = 0.2


time_Myr, sdf = read_dump('cloud_20_10_clrsink14_02900')

ihighrho = np.where(sdf['rho'] > 1e-23)[0]
sdf = sdf.iloc[ihighrho]

islice = np.where(np.logical_and(sdf['z'] > zmin1,sdf['z'] < zmax1))[0]
sdf_slice1 = sdf.iloc[islice]
plot_render_u(ax1,time_Myr,sdf_slice1,xmin,xmax,ymin,ymax,zmin1,zmax1,u_min,u_max)

islice = np.where(np.logical_and(sdf['z'] > zmin2,sdf['z'] < zmax2))[0]
sdf_slice2 = sdf.iloc[islice]
plot_render_u(ax2,time_Myr,sdf_slice2,xmin,xmax,ymin,ymax,zmin2,zmax2,u_min,u_max)



n,i,x,y,z,h,m,nH = np.loadtxt('nixyzhmf_00010.txt',skiprows=3,unpack=True)

pcm = plot_slice_nH(ax3,time_Myr,xmin,xmax,ymin,ymax,zmin1,zmax1)
pcm = plot_slice_nH(ax4,time_Myr,xmin,xmax,ymin,ymax,zmin2,zmax2)

ax3.set_aspect('equal','box')
ax4.set_aspect('equal','box')



cbar = fig.colorbar(pcm, ax=ax3, location='right', shrink=0.9,pad=0.02)
cbar.set_label('neutral fraction',rotation=270,labelpad=15)

cbar = fig.colorbar(pcm, ax=ax4, location='right', shrink=0.9,pad=0.02)
cbar.set_label('neutral fraction',rotation=270,labelpad=15)

fig.tight_layout(pad=0.7)


plt.savefig('hiiregion.png',dpi=200)
plt.savefig('hiiregion.pdf',format='pdf')
plt.show() 

