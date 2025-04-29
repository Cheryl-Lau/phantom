
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.cm import ScalarMappable
import sarracen

import warnings
warnings.filterwarnings('ignore')

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

    sdf_sinks = sdf_sinks[sdf_sinks['m'] > 0]

    return time_Myr, sdf, sdf_sinks



def get_neigh_stellarmass(sdf_sinks):

    rad_thresh = 1.0  # pc

    sdf_sinks['neigh_mass'] = np.zeros(len(sdf_sinks.index))

    for isink, sdf_isink in sdf_sinks.iterrows():
        pos_isink = np.array([sdf_isink['x'],sdf_isink['y'],sdf_isink['z']])
        for jsink, sdf_jsink in sdf_sinks.iterrows():
            pos_jsink = np.array([sdf_jsink['x'],sdf_jsink['y'],sdf_jsink['z']])
            m_jsink = sdf_jsink['m']
            if (isink != jsink):
                dist = np.linalg.norm(pos_isink - pos_jsink)
                if (dist < rad_thresh):
                    sdf_sinks['neigh_mass'][isink] += m_jsink
        
    print(sdf_sinks['neigh_mass'])

    return sdf_sinks



def plot_imf(ax,sdf_sinks,rho0str,alpha,time_Myr,t_ff):
    
    mass_min = 1e-1
    mass_max = 1e+2

    N = 50
    logmass_min = np.log10(mass_min)
    logmass_max = np.log10(mass_max)
    m_binned = np.logspace(logmass_min,logmass_max,N)
    dm = (logmass_max-logmass_min)/N 

    # Add extra column in sdf_sinks to count neighbouring stellar mass 
    sdf_sinks = get_neigh_stellarmass(sdf_sinks)
    
    N_binned = [0]*N
    totmassneigh_binned = [0.]*N

    for isink, sdf_isink in sdf_sinks.iterrows():
        m = sdf_isink['m'] 
        imbin = int((np.log10(m) - logmass_min + dm) /dm)-1
        if (imbin >= N):
            imbin = N-1
        elif (imbin < 0):
            imbin = 0
        N_binned[imbin] += 1 
        totmassneigh_binned[imbin] += sdf_isink['neigh_mass']
        
    # Get average neigh-mass in each bin
    totmassneigh_binned = np.divide(np.array(totmassneigh_binned), np.array(N_binned), out=np.zeros_like(np.array(totmassneigh_binned)), where=np.array(N_binned)!=0)    

    print(totmassneigh_binned)

    # Normalize to max 1 
    maxN = max(N_binned)
    for i in range(len(N_binned)):
        if (N_binned[i] > 0):
            N_binned[i] = N_binned[i]/maxN

    # Create colourbar 
    totm_color_norm = [x / max(totmassneigh_binned) for x in totmassneigh_binned]
    my_cmap = plt.cm.get_cmap('Blues')
    colors = my_cmap(totm_color_norm)

    ax.bar(m_binned[:-1], np.array(N_binned[:-1]),width=np.diff(m_binned),log=True,align='edge',color=colors)

    sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(min(totmassneigh_binned),max(totmassneigh_binned)))
    cbar = plt.colorbar(sm,ax=ax)
    cbar.set_label('Mean neighbouring sink mass [$M_\odot$]', rotation=270,labelpad=12)

    plot_settings(ax)

    return ax


def plot_mass_distri(ax,sdf_sinks):

    mass_min = 1e-1
    mass_max = 1e+2

    N = 50
    masses = sdf_sinks['m']
    ax.hist(masses,np.logspace(np.log10(mass_min),np.log10(mass_max), N))

    ax.set_xscale('log')
    ax.set_yscale('log')

    return ax


def plot_settings(ax):

    ax.set_xlabel('sink mass [$M_\odot$]')
    ax.set_ylabel('d$N$/dlog($m$)')
    ax.set_xscale('log')
    ax.set_yscale('log')

    xmin = 7e-2
    xmax = 2e+2
    ymin = 1e-2 
    ymax = 1.8e+0
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.text(1e-1,1.2e+0,'t = '+str(round(time_Myr,3))+' Myr ('+str(round(time_Myr/t_ff,3))+' $t_{ff}$)',color='black')
    ax.text(1e-1,1.5e-2,rho0str+r'$\ \ \ \alpha =$'+str(format(alpha,'.2f')),color='black',bbox=dict(facecolor='white',alpha=0.5))

    return 


fig, axes = plt.subplots(nrows=3, ncols=3, subplot_kw=dict(frameon=True), figsize=(12,10))


##### rho 1e-22 #####

rho0str = r'$\rho_0 = 10^{-22}\ \mathrm{g\ cm^{-3}}$'
t_ff = 6.662 # Myr 

ax22_07 = axes[0,0]
alpha = 0.7
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_07/cloud_22_07_07025')
plot_imf(ax22_07,sdf_sinks,rho0str,alpha,time_Myr,t_ff)

ax22_10 = axes[1,0]
alpha = 1.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_10/cloud_22_10_08100')
plot_imf(ax22_10,sdf_sinks,rho0str,alpha,time_Myr,t_ff)

ax22_20 = axes[2,0]
alpha = 2.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_20/cloud_22_20_08700')
plot_imf(ax22_20,sdf_sinks,rho0str,alpha,time_Myr,t_ff)


##### rho 1e-21 #####

rho0str = r'$\rho_0 = 10^{-21}\ \mathrm{g\ cm^{-3}}$'
t_ff = 2.107 # Myr 

ax21_07 = axes[0,1]
alpha = 0.7
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_07/cloud_21_07_02462')
plot_imf(ax21_07,sdf_sinks,rho0str,alpha,time_Myr,t_ff)

ax21_10 = axes[1,1]
alpha = 1.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_10/cloud_21_10_12010')
plot_imf(ax21_10,sdf_sinks,rho0str,alpha,time_Myr,t_ff)

ax21_20 = axes[2,1]
alpha = 2.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_20/cloud_21_20_10900')
plot_imf(ax21_20,sdf_sinks,rho0str,alpha,time_Myr,t_ff)


##### rho 1e-20 #####

rho0str = r'$\rho_0 = 10^{-20}\ \mathrm{g\ cm^{-3}}$'
t_ff = 0.666 # Myr 

ax20_07 = axes[0,2]
alpha = 0.7
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_07/cloud_20_07_00585')
plot_imf(ax20_07,sdf_sinks,rho0str,alpha,time_Myr,t_ff)

ax20_10 = axes[1,2]
alpha = 1.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_10/cloud_20_10_03748')
plot_imf(ax20_10,sdf_sinks,rho0str,alpha,time_Myr,t_ff)

ax20_20 = axes[2,2]
alpha = 2.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_20/cloud_20_20_02800')
plot_imf(ax20_20,sdf_sinks,rho0str,alpha,time_Myr,t_ff)


fig.tight_layout()
fig.subplots_adjust(wspace=0.32)
fig.subplots_adjust(hspace=0.32)

plt.savefig('sink_mass_distribution.png',dpi=200)
plt.savefig('sink_mass_distribution.pdf',format='pdf')



'''

fig2, axes2 = plt.subplots(nrows=3, ncols=3, subplot_kw=dict(frameon=True), figsize=(12,12))

##### rho 1e-22 #####

rho0str = r'$\rho_0 = 10^{-22}\ \mathrm{g\ cm^{-3}}$'

ax22_07 = axes2[0,0]
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_07/cloud_22_07_07025')
plot_mass_distri(ax22_07,sdf_sinks)

ax22_10 = axes2[1,0]
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_10/cloud_22_10_08100')
plot_mass_distri(ax22_10,sdf_sinks)

ax22_20 = axes2[2,0]
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_20/cloud_22_20_08700')
plot_mass_distri(ax22_20,sdf_sinks)


##### rho 1e-21 #####

rho0str = r'$\rho_0 = 10^{-21}\ \mathrm{g\ cm^{-3}}$'

ax21_07 = axes2[0,1]
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_07/cloud_21_07_02462')
plot_mass_distri(ax21_07,sdf_sinks)

ax21_10 = axes2[1,1]
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_10/cloud_21_10_12010')
plot_mass_distri(ax21_10,sdf_sinks)

ax21_20 = axes2[2,1]
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_20/cloud_21_20_10900')
plot_mass_distri(ax21_20,sdf_sinks)


##### rho 1e-20 #####

rho0str = r'$\rho_0 = 10^{-20}\ \mathrm{g\ cm^{-3}}$'

ax20_07 = axes2[0,2]
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_07/cloud_20_07_00585')
plot_mass_distri(ax20_07,sdf_sinks)

ax20_10 = axes2[1,2]
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_10/cloud_20_10_03748')
plot_mass_distri(ax20_10,sdf_sinks)

ax20_20 = axes2[2,2]
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_20/cloud_20_20_02800')
plot_mass_distri(ax20_20,sdf_sinks)

'''


plt.show()


























