

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

kboltz = 1.38066e-16
gmw = 2.381 
mass_proton_cgs = 1.67262158e-24


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
    sdf['T'] = sdf['u']/kboltz*(gmw*mass_proton_cgs*(gamma-1.))

    return time_Myr, sdf, sdf_sinks


def plot_rho_vs_radius(ax,sdf,x_src,y_src,z_src):

    x = sdf['x']
    y = sdf['y']
    z = sdf['z']
    r = np.array([x,y,z]).T

    r_src0 = np.array([x_src,y_src,z_src])
    r_src = np.stack((r_src0,) * len(sdf), axis=0)

    dr = r - r_src   
    absdr = np.linalg.norm(dr, axis=1)
    sdf['dr'] = absdr 

    as_colourbar = True

    if (as_colourbar):

        cm = plt.cm.get_cmap('inferno')
        pcm = ax.scatter(sdf['dr'],sdf['rho'],s=1,c=np.log10(sdf['T']),alpha=0.4,vmin=0.8,vmax=5.2,cmap=cm)

        cbar = fig.colorbar(pcm, ax=ax, location='right', shrink=0.8,pad=0.06)
        cbar.set_label('log T [K]',rotation=270,labelpad=10)
    else: 

        sdf_hot = sdf[sdf['T'] > 10000]
        x_hot = sdf_hot['x']
        y_hot = sdf_hot['y']
        z_hot = sdf_hot['z']
        r_hot = np.array([x_hot,y_hot,z_hot]).T

        r_src = np.stack((r_src0,) * len(sdf_hot), axis=0)

        dr_hot = r_hot - r_src   
        absdr_hot = np.linalg.norm(dr_hot, axis=1)
        sdf_hot['dr'] = absdr_hot

        ax.scatter(sdf_hot['dr'],sdf_hot['rho'],s=1,c='gold',alpha=1,label='T > 10000 K')
        

        sdf_cold = sdf[sdf['T'] < 50]
        x_cold = sdf_cold['x']
        y_cold = sdf_cold['y']
        z_cold = sdf_cold['z']
        r_cold = np.array([x_cold,y_cold,z_cold]).T

        r_src = np.stack((r_src0,) * len(sdf_cold), axis=0)

        dr_cold = r_cold - r_src   
        absdr_cold = np.linalg.norm(dr_cold, axis=1)
        sdf_cold['dr'] = absdr_cold

        ax.scatter(sdf_cold['dr'],sdf_cold['rho'],s=1,c='navy',alpha=1,label='T < 50 K')

    return 


def get_most_massive_sink(sdf_sinks):

    m_max = np.max(sdf_sinks['m'])
    imax = sdf_sinks.index[sdf_sinks['m'] == m_max] 
    
    isink_max = imax[0] + 1 
    x_max = sdf_sinks['x'].loc[imax]
    y_max = sdf_sinks['y'].loc[imax]
    z_max = sdf_sinks['z'].loc[imax]

    return isink_max, x_max.tolist()[0], y_max.tolist()[0], z_max.tolist()[0]



def plot_settings(ax,time_Myr,t_ion,rho0str,alpha):

    xmin = -0.6
    xmax = 7.1
    ymin = 1e-25
    ymax = 1e-14

    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])

    ax.set_yscale('log')

    ax.set_xlabel('r [pc]')
    ax.set_ylabel('density '+r'$\mathrm{[g\ cm^{-3}]}$')
    #ax.text(xmin + 0.6*(xmax-xmin), ymax*0.001,rho0str,color='black',fontsize=8)
    #ax.text(xmin + 0.8*(xmax-xmin), ymax*0.0005,r'$\alpha = $'+str(format(alpha,'.2f')),fontsize=8)

    ax.legend(loc='upper right',fontsize=8)

    return 


def plot_cavline(ax,r_cav):

    ax.axvline(r_cav, linestyle='--', linewidth=0.5, color='royalblue', label='cavity')

    return 



'''
With channels 
'''

fig, axes = plt.subplots(nrows=2, ncols=2, subplot_kw=dict(frameon=True), figsize=(8,6))


##### rho 1e-21 #####

rho0str = r'$\rho_0 = 10^{-21}\ \mathrm{g\ cm^{-3}}$'
t_ff = 2.107 # Myr 

ax21_07 = axes[0,0]
alpha = 0.7
time0_Myr, sdf, sdf_sinks = read_dump('rad_21_07/cloud_21_07_clrsink4_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_21_07/cloud_21_07_clrsink4_28572')  
t_ion = time_Myr - time0_Myr 
plot_cavline(ax21_07,0.36)
plot_rho_vs_radius(ax21_07,sdf,x_src,y_src,z_src)
plot_settings(ax21_07,time_Myr,t_ion,rho0str,alpha)

ax21_10 = axes[1,0]
alpha = 1.0
time0_Myr, sdf, sdf_sinks = read_dump('rad_21_10/cloud_21_10_clrsink6_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_21_10/cloud_21_10_clrsink6_25410')  
t_ion = time_Myr - time0_Myr 
plot_cavline(ax21_10,0.53)
plot_rho_vs_radius(ax21_10,sdf,x_src,y_src,z_src)
plot_settings(ax21_10,time_Myr,t_ion,rho0str,alpha)



##### rho 1e-20 #####

rho0str = r'$\rho_0 = 10^{-20}\ \mathrm{g\ cm^{-3}}$'
t_ff = 0.666 # Myr 

ax20_07 = axes[0,1]
alpha = 0.7
time0_Myr, sdf, sdf_sinks = read_dump('rad_20_07/cloud_20_07_clrsink1_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_20_07/cloud_20_07_clrsink1_21490') 
t_ion = time_Myr - time0_Myr 
plot_cavline(ax20_07,0.238)
plot_rho_vs_radius(ax20_07,sdf,x_src,y_src,z_src)
plot_settings(ax20_07,time_Myr,t_ion,rho0str,alpha)

ax20_10 = axes[1,1]
alpha = 1.0
time0_Myr, sdf, sdf_sinks = read_dump('rad_20_10/cloud_20_10_clrsink14_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_20_10/cloud_20_10_clrsink14_20472') 
t_ion = time_Myr - time0_Myr 
plot_cavline(ax20_10,0.268)
plot_rho_vs_radius(ax20_10,sdf,x_src,y_src,z_src)
plot_settings(ax20_10,time_Myr,t_ion,rho0str,alpha)



fig.tight_layout()
fig.subplots_adjust(hspace=0.03)

plt.savefig('chnl_density_intermed.png',dpi=200)



'''
Without channels 
'''

fig, axes = plt.subplots(nrows=2, ncols=2, subplot_kw=dict(frameon=True), figsize=(8,6))


##### rho 1e-21 #####

rho0str = r'$\rho_0 = 10^{-21}\ \mathrm{g\ cm^{-3}}$'
t_ff = 2.107 # Myr 

ax21_07 = axes[0,0]
alpha = 0.7
time0_Myr, sdf, sdf_sinks = read_dump('rad_21_07/cloud_21_07_clrsink4_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_21_07/cloud_21_07_clrsink4_29329')  
t_ion = time_Myr - time0_Myr 
plot_cavline(ax21_07,2.97)
plot_rho_vs_radius(ax21_07,sdf,x_src,y_src,z_src)
plot_settings(ax21_07,time_Myr,t_ion,rho0str,alpha)

ax21_10 = axes[1,0]
alpha = 1.0
time0_Myr, sdf, sdf_sinks = read_dump('rad_21_10/cloud_21_10_clrsink6_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_21_10/cloud_21_10_clrsink6_26245')  
t_ion = time_Myr - time0_Myr 
plot_cavline(ax21_10,4.18)
plot_rho_vs_radius(ax21_10,sdf,x_src,y_src,z_src)
plot_settings(ax21_10,time_Myr,t_ion,rho0str,alpha)


##### rho 1e-20 #####

rho0str = r'$\rho_0 = 10^{-20}\ \mathrm{g\ cm^{-3}}$'
t_ff = 0.666 # Myr 

ax20_07 = axes[0,1]
alpha = 0.7
time0_Myr, sdf, sdf_sinks = read_dump('rad_20_07/cloud_20_07_clrsink1_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_20_07/cloud_20_07_clrsink1_24172') 
t_ion = time_Myr - time0_Myr 
plot_cavline(ax20_07,1.54)
plot_rho_vs_radius(ax20_07,sdf,x_src,y_src,z_src)
plot_settings(ax20_07,time_Myr,t_ion,rho0str,alpha)

ax20_10 = axes[1,1]
alpha = 1.0
time0_Myr, sdf, sdf_sinks = read_dump('rad_20_10/cloud_20_10_clrsink14_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_20_10/cloud_20_10_clrsink14_22955') 
t_ion = time_Myr - time0_Myr 
plot_cavline(ax20_10,1.44)
plot_rho_vs_radius(ax20_10,sdf,x_src,y_src,z_src)
plot_settings(ax20_10,time_Myr,t_ion,rho0str,alpha)



fig.tight_layout()
fig.subplots_adjust(hspace=0.03)

plt.savefig('chnl_density_final.png',dpi=200)


plt.show()