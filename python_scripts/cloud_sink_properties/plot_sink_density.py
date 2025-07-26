
import numpy as np
import matplotlib.pyplot as plt 
import sarracen
from scipy.optimize import curve_fit

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



def get_neigh_stellarmass(sdf_sinks,radin_thresh,radout_thresh):

    sdf_sinks['neigh_mass'] = np.zeros(len(sdf_sinks.index))
    sdf_sinks['nneigh'] = np.zeros(len(sdf_sinks.index))

    for isink, sdf_isink in sdf_sinks.iterrows():
        pos_isink = np.array([sdf_isink['x'],sdf_isink['y'],sdf_isink['z']])
        for jsink, sdf_jsink in sdf_sinks.iterrows():
            pos_jsink = np.array([sdf_jsink['x'],sdf_jsink['y'],sdf_jsink['z']])
            m_jsink = sdf_jsink['m']
            if (isink != jsink):
                dist = np.linalg.norm(pos_isink - pos_jsink)
                if (dist > radin_thresh and dist < radout_thresh):
                    sdf_sinks['nneigh'][isink] += 1
                    sdf_sinks['neigh_mass'][isink] += m_jsink
        
    # Calculate mass density and number density 
    vol_thresh = 4.0/3.0*np.pi*radout_thresh**3 - 4.0/3.0*np.pi*radin_thresh**3
    sdf_sinks['rho'] = sdf_sinks['neigh_mass']/vol_thresh
    sdf_sinks['nrho'] = sdf_sinks['nneigh']/vol_thresh

    # Calculate mean neighbour mass 
    sdf_sinks['avg_neigh_mass'] = sdf_sinks['neigh_mass']/sdf_sinks['nneigh']

    return sdf_sinks



def sink_loc_in_cloud(sdf_sinks,NS_divide):
    
    sdf_sinks['loc_in_cloud'] = ['O']*len(sdf_sinks.index)

    iNorth = np.where(sdf_sinks['y'] > NS_divide)[0]
    sdf_sinks['loc_in_cloud'].iloc[iNorth] = 'N'
    
    iSouth = np.where(sdf_sinks['y'] <= NS_divide)[0]
    sdf_sinks['loc_in_cloud'].iloc[iSouth] = 'S'

    return sdf_sinks



def plot_sink_mass_neighmass(ax,sdf_sinks,radin_thresh):

    radout_thresh = radin_thresh+0.5  # pc
    sdf_sinks = get_neigh_stellarmass(sdf_sinks,radin_thresh,radout_thresh)
    sdf_sinks_N = sdf_sinks[sdf_sinks['loc_in_cloud'] == 'N']
    sdf_sinks_S = sdf_sinks[sdf_sinks['loc_in_cloud'] == 'S']
    ax.scatter(sdf_sinks_N['m'],sdf_sinks_N['rho'],s=2,color='black',label=str(round(radin_thresh,2))+' < r < '+str(round(radout_thresh,2))+' pc in N')
    ax.scatter(sdf_sinks_S['m'],sdf_sinks_S['rho'],s=10,facecolors='none',edgecolors='black',label=str(round(radin_thresh,2))+' < r < '+str(round(radout_thresh,2))+' pc in S')

    radout_thresh = radin_thresh+1.0  # pc
    sdf_sinks = get_neigh_stellarmass(sdf_sinks,radin_thresh,radout_thresh)
    sdf_sinks_N = sdf_sinks[sdf_sinks['loc_in_cloud'] == 'N']
    sdf_sinks_S = sdf_sinks[sdf_sinks['loc_in_cloud'] == 'S']
    ax.scatter(sdf_sinks_N['m'],sdf_sinks_N['rho'],s=2,color='blue',label=str(round(radin_thresh,2))+' < r < '+str(round(radout_thresh,2))+' pc in N')
    ax.scatter(sdf_sinks_S['m'],sdf_sinks_S['rho'],s=10,facecolors='none',edgecolors='blue',label=str(round(radin_thresh,2))+' < r < '+str(round(radout_thresh,2))+' pc in S')

    radout_thresh = radin_thresh+1.5  # pc
    sdf_sinks = get_neigh_stellarmass(sdf_sinks,radin_thresh,radout_thresh)
    sdf_sinks_N = sdf_sinks[sdf_sinks['loc_in_cloud'] == 'N']
    sdf_sinks_S = sdf_sinks[sdf_sinks['loc_in_cloud'] == 'S']
    ax.scatter(sdf_sinks_N['m'],sdf_sinks_N['rho'],s=2,color='red',label=str(round(radin_thresh,2))+' < r < '+str(round(radout_thresh,2))+' pc in N')
    ax.scatter(sdf_sinks_S['m'],sdf_sinks_S['rho'],s=10,facecolors='none',edgecolors='red',label=str(round(radin_thresh,2))+' < r < '+str(round(radout_thresh,2))+' pc in S')

    radout_thresh = radin_thresh+2.0  # pc
    sdf_sinks = get_neigh_stellarmass(sdf_sinks,radin_thresh,radout_thresh)
    sdf_sinks_N = sdf_sinks[sdf_sinks['loc_in_cloud'] == 'N']
    sdf_sinks_S = sdf_sinks[sdf_sinks['loc_in_cloud'] == 'S']
    ax.scatter(sdf_sinks_N['m'],sdf_sinks_N['rho'],s=2,color='green',label=str(round(radin_thresh,2))+' < r < '+str(round(radout_thresh,2))+' pc in N')
    ax.scatter(sdf_sinks_S['m'],sdf_sinks_S['rho'],s=10,facecolors='none',edgecolors='green',label=str(round(radin_thresh,2))+' < r < '+str(round(radout_thresh,2))+' pc in S')

    return 


def plot_fit_line(ax,x,y):

    param, covar = curve_fit(linear_line, x, y, p0=[-1,200])
    k = param[0]
    C = param[1]

    mass = np.linspace(0,160,1000)
    ax.plot(mass,linear_line(mass, -1, C),'--',color='grey',label=r'$m_\mathrm{neigh} + m = C$')

    return 


def plot_line(ax,a):

    mass = np.logspace(0,3,1000)
    ax.plot(mass,power_law(mass, a, -1),'--',color='lightgrey',label=r'$\rho_\mathrm{neigh} \propto m^{-1}$')

    return 


def linear_line(x, k, C):

    y = k*x + C

    return y

def power_law(x, a, k):

    y = a*x**k

    return y


def plot_settings(ax,time_Myr,t_ff,rho0str,alpha):

    ax.set_xlabel('sink mass $m$ [$M_\odot$]')
    ax.set_ylabel(r'$\rho_\mathrm{neigh} $ [$M_\odot$ $\mathrm{pc}^{-3}$]')   

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_xlim([8e0,4e2])
    ax.set_ylim([1.9e-2,2e2])
#    ax.text(xmin+(xmax-xmin)*0.03,ymax-(ymax-ymin)*0.1,rho0str,color='black',fontsize=9)
#    ax.text(xmin+(xmax-xmin)*0.03,ymax-(ymax-ymin)*0.2,r'$\alpha =$'+str(format(alpha,'.2f')),color='black',fontsize=9)

    ax.legend(fontsize=5,loc='upper right')

    return 



fig, axes = plt.subplots(nrows=3, ncols=3, subplot_kw=dict(frameon=True), figsize=(10,9))


##### rho 1e-22 #####

rho0str = r'$\rho_0 = 10^{-22}\ \mathrm{g\ cm^{-3}}$'
t_ff = 6.662 # Myr 
NS_divide = 7.0

radin_thresh = 1.3

ax22_07 = axes[0,0]
alpha = 0.7
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_07/cloud_22_07_07025')
sdf_sinks = sink_loc_in_cloud(sdf_sinks,NS_divide)
plot_sink_mass_neighmass(ax22_07,sdf_sinks,radin_thresh)
plot_line(ax22_07,2e2)
plot_settings(ax22_07,time_Myr,t_ff,rho0str,alpha)

radin_thresh = 1.2

ax22_10 = axes[1,0]
alpha = 1.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_10/cloud_22_10_08100')
sdf_sinks = sink_loc_in_cloud(sdf_sinks,NS_divide)
plot_sink_mass_neighmass(ax22_10,sdf_sinks,radin_thresh)
plot_line(ax22_10,5e2)
plot_settings(ax22_10,time_Myr,t_ff,rho0str,alpha)

radin_thresh = 1.5

ax22_20 = axes[2,0]
alpha = 2.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_22_20/cloud_22_20_08700')
sdf_sinks = sink_loc_in_cloud(sdf_sinks,NS_divide)
plot_sink_mass_neighmass(ax22_20,sdf_sinks,radin_thresh)
plot_line(ax22_20,2e1)
plot_settings(ax22_20,time_Myr,t_ff,rho0str,alpha)



##### rho 1e-21 #####

rho0str = r'$\rho_0 = 10^{-21}\ \mathrm{g\ cm^{-3}}$'
t_ff = 2.107 # Myr 
NS_divide = 3.0

radin_thresh = 1.6  

ax21_07 = axes[0,1]
alpha = 0.7
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_07/cloud_21_07_02462')
sdf_sinks = sink_loc_in_cloud(sdf_sinks,NS_divide)
plot_sink_mass_neighmass(ax21_07,sdf_sinks,radin_thresh)
plot_line(ax21_07,2e2)
plot_settings(ax21_07,time_Myr,t_ff,rho0str,alpha)

radin_thresh = 1.7

ax21_10 = axes[1,1]
alpha = 1.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_10/cloud_21_10_12010')
sdf_sinks = sink_loc_in_cloud(sdf_sinks,NS_divide)
plot_sink_mass_neighmass(ax21_10,sdf_sinks,radin_thresh)
plot_line(ax21_10,2e2)
plot_settings(ax21_10,time_Myr,t_ff,rho0str,alpha)

radin_thresh = 3.5

ax21_20 = axes[2,1]
alpha = 2.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_21_20/cloud_21_20_10900')
sdf_sinks = sink_loc_in_cloud(sdf_sinks,NS_divide)
plot_sink_mass_neighmass(ax21_20,sdf_sinks,radin_thresh)
plot_line(ax21_20,2e2)
plot_settings(ax21_20,time_Myr,t_ff,rho0str,alpha)



##### rho 1e-20 #####

rho0str = r'$\rho_0 = 10^{-20}\ \mathrm{g\ cm^{-3}}$'
t_ff = 0.666 # Myr 
NS_divide = 1.0

radin_thresh = 0.3

ax20_07 = axes[0,2]
alpha = 0.7
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_07/cloud_20_07_00585')
sdf_sinks = sink_loc_in_cloud(sdf_sinks,NS_divide)
plot_sink_mass_neighmass(ax20_07,sdf_sinks,radin_thresh)
plot_line(ax20_07,2e2)
plot_settings(ax20_07,time_Myr,t_ff,rho0str,alpha)

radin_thresh = 0.5

ax20_10 = axes[1,2]
alpha = 1.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_10/cloud_20_10_03748')
sdf_sinks = sink_loc_in_cloud(sdf_sinks,NS_divide)
plot_sink_mass_neighmass(ax20_10,sdf_sinks,radin_thresh)
plot_line(ax20_10,4e2)
plot_settings(ax20_10,time_Myr,t_ff,rho0str,alpha)

radin_thresh = 0.5

ax20_20 = axes[2,2]
alpha = 2.0
time_Myr, sdf, sdf_sinks = read_dump('elongev_20_20/cloud_20_20_02800')
sdf_sinks = sink_loc_in_cloud(sdf_sinks,NS_divide)
plot_sink_mass_neighmass(ax20_20,sdf_sinks,radin_thresh)
plot_line(ax20_20,4e2)
plot_settings(ax20_20,time_Myr,t_ff,rho0str,alpha)



fig.tight_layout()


plt.savefig('sink_neighrho.png',dpi=200)
plt.savefig('sink_neighrho.pdf',format='pdf')
plt.show()










