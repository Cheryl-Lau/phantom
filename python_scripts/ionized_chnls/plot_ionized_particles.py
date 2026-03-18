
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

    return time_Myr, sdf, sdf_sinks


def plot_slice_nH(ax,x_source,y_source,ionfrac_outfile,zmin,zmax):

    n,i,x,y,z,h,m,nH = np.loadtxt(ionfrac_outfile,skiprows=3,unpack=True)

    islice = np.where((z > zmin) & (z < zmax))[0]
    x_slice = x[islice]
    y_slice = y[islice]
    m_slice = m[islice]
    h_slice = h[islice]
    nH_slice = nH[islice]

    pcm = ax.scatter(x_slice,y_slice,s=1,c=nH_slice,alpha=0.4)
    
    cbar = fig.colorbar(pcm, ax=ax, location='right', shrink=0.8,pad=0.06)
    cbar.set_label('neutral fraction',rotation=270,labelpad=10)

    # Mark the ionizing sink  
    ax.scatter(x_source,y_source,s=70,marker='*',color='salmon',label='ionizing source')
    ax.legend(loc='lower right')

    return 


def plot_settings(ax,time_Myr,t_ion,rho0str,alpha,p_chnl,xmin,xmax,ymin,ymax):

    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.9*(ymax-ymin),rho0str,color='black')
    ax.text(xmax - 0.4*(xmax-xmin), ymin + 0.9*(ymax-ymin),r'$\alpha = $'+str(format(alpha,'.2f')))

    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.83*(ymax-ymin),r'$t_\mathrm{cloud} = $'+str(format(time_Myr,'.2f'))+' Myr')
    ax.text(xmax - 0.4*(xmax-xmin), ymin + 0.83*(ymax-ymin),r'$t_\mathrm{ion} = $'+str(format(t_ion,'.3f'))+' Myr')

    ax.text(xmax - 0.52*(xmax-xmin), ymin + 0.15*(ymax-ymin),str(round(zmin,1))+' pc < z < '+str(round(zmax,1))+' pc')

    ax.text(xmin + 0.08*(xmax-xmin), ymin + 0.05*(ymax-ymin),r'$P_\mathrm{chnl} = $'+str(p_chnl))

    ax.set_aspect('equal','box')

    return 


def get_most_massive_sink(ax,sdf_sinks):

    m_max = np.max(sdf_sinks['m'])
    imax = sdf_sinks.index[sdf_sinks['m'] == m_max] 
    
    isink_max = imax[0] + 1 
    x_max = sdf_sinks['x'].loc[imax]
    y_max = sdf_sinks['y'].loc[imax]
    z_max = sdf_sinks['z'].loc[imax]

    return isink_max, x_max.tolist()[0], y_max.tolist()[0], z_max.tolist()[0]


'''
With channels 
'''

fig, axes = plt.subplots(nrows=2, ncols=2, subplot_kw=dict(frameon=True), figsize=(9,8))


##### rho 1e-21 #####

rho0str = r'$\rho_0 = 10^{-21}\ \mathrm{g\ cm^{-3}}$'
t_ff = 2.107 # Myr 

ax21_07 = axes[0,0]
alpha = 0.7
p_chnl = 128.5
time0_Myr, sdf, sdf_sinks = read_dump('rad_21_07/cloud_21_07_clrsink4_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(ax21_07,sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_21_07/cloud_21_07_clrsink4_28572')  
t_ion = time_Myr - time0_Myr 
radius = 8.0
xmin = x_src - radius 
xmax = x_src + radius 
ymin = y_src - radius 
ymax = y_src + radius 
zmin = z_src - 0.05
zmax = z_src + 0.05
print('isink 21_07: ',isink)
plot_slice_nH(ax21_07,x_src,y_src,'rad_21_07/new_nixyzhmf/nixyzhmf_00003.txt',zmin,zmax) 
plot_settings(ax21_07,time_Myr,t_ion,rho0str,alpha,p_chnl,xmin,xmax,ymin,ymax)

ax21_10 = axes[1,0]
alpha = 1.0
p_chnl = 82.8
time0_Myr, sdf, sdf_sinks = read_dump('rad_21_10/cloud_21_10_clrsink6_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(ax21_10,sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_21_10/cloud_21_10_clrsink6_25410')  
t_ion = time_Myr - time0_Myr 
radius = 8.0
xmin = x_src - radius 
xmax = x_src + radius 
ymin = y_src - radius 
ymax = y_src + radius 
zmin = z_src - 0.05
zmax = z_src + 0.05
print('isink 21_10: ',isink)
plot_slice_nH(ax21_10,x_src,y_src,'rad_21_10/new_nixyzhmf/nixyzhmf_00003.txt',zmin,zmax) 
plot_settings(ax21_10,time_Myr,t_ion,rho0str,alpha,p_chnl,xmin,xmax,ymin,ymax)


##### rho 1e-20 #####

rho0str = r'$\rho_0 = 10^{-20}\ \mathrm{g\ cm^{-3}}$'
t_ff = 0.666 # Myr 

ax20_07 = axes[0,1]
alpha = 0.7
p_chnl = 306.5
time0_Myr, sdf, sdf_sinks = read_dump('rad_20_07/cloud_20_07_clrsink1_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(ax20_07,sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_20_07/cloud_20_07_clrsink1_21490') 
t_ion = time_Myr - time0_Myr 
radius = 8.0
xmin = x_src - radius 
xmax = x_src + radius 
ymin = y_src - radius 
ymax = y_src + radius 
zmin = z_src - 0.05
zmax = z_src + 0.05
print('isink 20_07: ',isink)
plot_slice_nH(ax20_07,x_src,y_src,'rad_20_07/new_nixyzhmf/nixyzhmf_00005.txt',zmin,zmax) 
plot_settings(ax20_07,time_Myr,t_ion,rho0str,alpha,p_chnl,xmin,xmax,ymin,ymax)

ax20_10 = axes[1,1]
alpha = 1.0
p_chnl = 174.2 
time0_Myr, sdf, sdf_sinks = read_dump('rad_20_10/cloud_20_10_clrsink14_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(ax20_10,sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_20_10/cloud_20_10_clrsink14_20472') 
t_ion = time_Myr - time0_Myr 
radius = 8.0
xmin = x_src - radius 
xmax = x_src + radius 
ymin = y_src - radius 
ymax = y_src + radius 
zmin = z_src - 0.05
zmax = z_src + 0.05
print('isink 20_10: ',isink)
plot_slice_nH(ax20_10,x_src,y_src,'rad_20_10/new_nixyzhmf/nixyzhmf_00020.txt',zmin,zmax) 
plot_settings(ax20_10,time_Myr,t_ion,rho0str,alpha,p_chnl,xmin,xmax,ymin,ymax)



fig.tight_layout()
fig.subplots_adjust(hspace=0.03)

plt.savefig('clouds_ionized_regions_intermed.png',dpi=200)



'''
Without channels 
'''

fig, axes = plt.subplots(nrows=2, ncols=2, subplot_kw=dict(frameon=True), figsize=(9,8))


##### rho 1e-21 #####

rho0str = r'$\rho_0 = 10^{-21}\ \mathrm{g\ cm^{-3}}$'
t_ff = 2.107 # Myr 

ax21_07 = axes[0,0]
alpha = 0.7
p_chnl = 4.0
time0_Myr, sdf, sdf_sinks = read_dump('rad_21_07/cloud_21_07_clrsink4_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(ax21_07,sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_21_07/cloud_21_07_clrsink4_29329')  # _28530
t_ion = time_Myr - time0_Myr 
radius = 8.0
xmin = x_src - radius 
xmax = x_src + radius 
ymin = y_src - radius 
ymax = y_src + radius 
zmin = z_src - 0.05
zmax = z_src + 0.05
print('isink 21_07: ',isink)
plot_slice_nH(ax21_07,x_src,y_src,'rad_21_07/new_nixyzhmf/nixyzhmf_00177.txt',zmin,zmax) # _00115
plot_settings(ax21_07,time_Myr,t_ion,rho0str,alpha,p_chnl,xmin,xmax,ymin,ymax)

ax21_10 = axes[1,0]
alpha = 1.0
p_chnl = 1.7
time0_Myr, sdf, sdf_sinks = read_dump('rad_21_10/cloud_21_10_clrsink6_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(ax21_10,sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_21_10/cloud_21_10_clrsink6_26245')  # _25350
t_ion = time_Myr - time0_Myr 
radius = 8.0
xmin = x_src - radius 
xmax = x_src + radius 
ymin = y_src - radius 
ymax = y_src + radius 
zmin = z_src - 0.05
zmax = z_src + 0.05
print('isink 21_10: ',isink)
plot_slice_nH(ax21_10,x_src,y_src,'rad_21_10/new_nixyzhmf/nixyzhmf_00098.txt',zmin,zmax) # _00255
plot_settings(ax21_10,time_Myr,t_ion,rho0str,alpha,p_chnl,xmin,xmax,ymin,ymax)


##### rho 1e-20 #####

rho0str = r'$\rho_0 = 10^{-20}\ \mathrm{g\ cm^{-3}}$'
t_ff = 0.666 # Myr 

ax20_07 = axes[0,1]
alpha = 0.7
p_chnl = 20.2
time0_Myr, sdf, sdf_sinks = read_dump('rad_20_07/cloud_20_07_clrsink1_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(ax20_07,sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_20_07/cloud_20_07_clrsink1_24172') # _21450
t_ion = time_Myr - time0_Myr 
radius = 8.0
xmin = x_src - radius 
xmax = x_src + radius 
ymin = y_src - radius 
ymax = y_src + radius 
zmin = z_src - 0.05
zmax = z_src + 0.05
print('isink 20_07: ',isink)
plot_slice_nH(ax20_07,x_src,y_src,'rad_20_07/new_nixyzhmf/nixyzhmf_00435.txt',zmin,zmax) # _00783
plot_settings(ax20_07,time_Myr,t_ion,rho0str,alpha,p_chnl,xmin,xmax,ymin,ymax)

ax20_10 = axes[1,1]
alpha = 1.0
p_chnl = 15.8
time0_Myr, sdf, sdf_sinks = read_dump('rad_20_10/cloud_20_10_clrsink14_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(ax20_10,sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_20_10/cloud_20_10_clrsink14_22955') # _20506
t_ion = time_Myr - time0_Myr 
radius = 8.0
xmin = x_src - radius 
xmax = x_src + radius 
ymin = y_src - radius 
ymax = y_src + radius 
zmin = z_src - 0.05
zmax = z_src + 0.05
print('isink 20_10: ',isink)
plot_slice_nH(ax20_10,x_src,y_src,'rad_20_10/new_nixyzhmf/nixyzhmf_00396.txt',zmin,zmax) # _00223
plot_settings(ax20_10,time_Myr,t_ion,rho0str,alpha,p_chnl,xmin,xmax,ymin,ymax)




fig.tight_layout()
fig.subplots_adjust(hspace=0.03)

plt.savefig('clouds_ionized_regions_final.png',dpi=200)
plt.show()




























