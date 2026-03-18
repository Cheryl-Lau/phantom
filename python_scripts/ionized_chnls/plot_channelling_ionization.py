

import numpy as np
import matplotlib.pyplot as plt 
import sarracen


utime = 4.706e14
udist = 3.086e+19
umass = 1.989e+33
unit_density = 6.768e-23
unit_velocity = 6.558e+3
unit_ergg = 4.301e+7
solarm = 1.989e+33
Myr = 1e6*365*24*60*60 

pmass_cgs = 1e-2*umass
pmass = pmass_cgs/solarm


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



def plot_disttosrc_hist(ax,x_src,y_src,z_src,ionfrac_outfile):

    n,i,x,y,z,h,m,nH = np.loadtxt(ionfrac_outfile,skiprows=3,unpack=True)

    r_src0 = np.array([x_src,y_src,z_src])

    # All particles 
    r = np.array([x,y,z]).T
    r_src = np.stack((r_src0,) * len(n), axis=0)
    dr = r - r_src   
    absdr = np.linalg.norm(dr, axis=1)

    counts, bins = np.histogram(absdr,200,range=(0,maxr))
    ax.stairs(counts*pmass, bins, color='red',label='total mass')


    # Extract neutral particles 
    ineu = np.where(nH > 0.5)[0]
    x_neu = x[ineu]
    y_neu = y[ineu]
    z_neu = z[ineu]
    r_neu = np.array([x_neu,y_neu,z_neu]).T

    r_src = np.stack((r_src0,) * len(ineu), axis=0)
    dr_neu = r_neu - r_src   
    absdr_neu = np.linalg.norm(dr_neu, axis=1)

    counts, bins = np.histogram(absdr_neu,200,range=(0,maxr))
    ax.stairs(counts*pmass, bins, color='gold',label='neutral mass')


    # Extract ionized particles 
    iion = np.where(nH < 0.5)[0]
    x_ion = x[iion]
    y_ion = y[iion]
    z_ion = z[iion]
    r_ion = np.array([x_ion,y_ion,z_ion]).T

    r_src = np.stack((r_src0,) * len(iion), axis=0)
    dr_ion = r_ion - r_src   
    absdr_ion = np.linalg.norm(dr_ion, axis=1)

    counts, bins = np.histogram(absdr_ion,200,range=(0,maxr))
    ax.stairs(counts*pmass, bins, color='purple',label='ionized mass')


    return 



def get_most_massive_sink(sdf_sinks):

    m_max = np.max(sdf_sinks['m'])
    imax = sdf_sinks.index[sdf_sinks['m'] == m_max] 
    
    isink_max = imax[0] + 1 
    x_max = sdf_sinks['x'].loc[imax]
    y_max = sdf_sinks['y'].loc[imax]
    z_max = sdf_sinks['z'].loc[imax]

    return isink_max, x_max.tolist()[0], y_max.tolist()[0], z_max.tolist()[0]


def plot_cavline(ax,r_cav):

    ax.axvline(r_cav, linestyle='--', linewidth=0.8, color='royalblue', label='cavity')

    return 


def plot_chnldist(ax,r_furthest):

    ax.axvline(r_furthest, linestyle='--', linewidth=0.8, color='lime', label='channel')

    return 



def plot_settings(ax,time_Myr,t_ion,rho0str,alpha):

    xmin = -0.5
    xmax = maxr*0.99
    ymin = 1.5e31/umass
    ymax = 8e35/umass

    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])

    ax.set_yscale('log')

    ax.set_xlabel('r [pc]')
    ax.set_ylabel('mass '+r'$\mathrm{[M_\odot]}$')
    ax.text(xmin + 0.6*(xmax-xmin), ymax*0.03,rho0str,color='black',fontsize=8)
    ax.text(xmin + 0.8*(xmax-xmin), ymax*0.015,r'$\alpha = $'+str(format(alpha,'.2f')),fontsize=8)

    ax.legend(loc='upper right',fontsize=6)

    return 



def get_ionizedmass_inshell(x_src,y_src,z_src,rmin,rmax,ionfrac_outfile):

    n,i,x,y,z,h,m,nH = np.loadtxt(ionfrac_outfile,skiprows=3,unpack=True)

    r_src = np.array([x_src,y_src,z_src])

    totmass = 0 
    for ip in range(len(n)): 
        nHi = nH[ip]
        if (nHi < 0.5):
            xi = x[ip]
            yi = y[ip]
            zi = z[ip]
            ri = np.array([xi,yi,zi])
            dri = ri - r_src
            absdri = np.linalg.norm(dri)
            if (absdri > rmin and absdri < rmax): 
                totmass += m[ip]

    return totmass


def get_neutralmass_inshell(x_src,y_src,z_src,rmin,rmax,ionfrac_outfile):

    n,i,x,y,z,h,m,nH = np.loadtxt(ionfrac_outfile,skiprows=3,unpack=True)

    r_src = np.array([x_src,y_src,z_src])

    totmass = 0 
    for ip in range(len(n)): 
        nHi = nH[ip]
        if (nHi > 0.5):
            xi = x[ip]
            yi = y[ip]
            zi = z[ip]
            ri = np.array([xi,yi,zi])
            dri = ri - r_src
            absdri = np.linalg.norm(dri)
            if (absdri > rmin and absdri < rmax): 
                totmass += m[ip]

    return totmass


'''
With channels 
'''

fig, axes = plt.subplots(nrows=2, ncols=2, subplot_kw=dict(frameon=True), figsize=(7,6))

maxr = 18.0 


##### rho 1e-21 #####

rho0str = r'$\rho_0 = 10^{-21}\ \mathrm{g\ cm^{-3}}$'
t_ff = 2.107 # Myr 

ax21_07 = axes[0,0]
alpha = 0.7
time0_Myr, sdf, sdf_sinks = read_dump('rad_21_07/cloud_21_07_clrsink4_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_21_07/cloud_21_07_clrsink4_28572')  
t_ion = time_Myr - time0_Myr 
plot_disttosrc_hist(ax21_07,x_src,y_src,z_src,'rad_21_07/new_nixyzhmf/nixyzhmf_00003.txt')
plot_cavline(ax21_07,0.36)
plot_chnldist(ax21_07,10.55)
plot_settings(ax21_07,time_Myr,t_ion,rho0str,alpha)
mcav = get_ionizedmass_inshell(x_src,y_src,z_src,0,0.36,'rad_21_07/new_nixyzhmf/nixyzhmf_00003.txt')
mchl_ion = get_ionizedmass_inshell(x_src,y_src,z_src,0.36,10.55,'rad_21_07/new_nixyzhmf/nixyzhmf_00003.txt')
mchl_neu = get_neutralmass_inshell(x_src,y_src,z_src,0.36,10.55,'rad_21_07/new_nixyzhmf/nixyzhmf_00003.txt')
print('21_07: ',mchl_ion,mchl_neu,mcav,mchl_neu/mcav)

ax21_10 = axes[1,0]
alpha = 1.0
time0_Myr, sdf, sdf_sinks = read_dump('rad_21_10/cloud_21_10_clrsink6_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_21_10/cloud_21_10_clrsink6_25410')  
t_ion = time_Myr - time0_Myr 
plot_disttosrc_hist(ax21_10,x_src,y_src,z_src,'rad_21_10/new_nixyzhmf/nixyzhmf_00003.txt') 
plot_cavline(ax21_10,0.53)
plot_chnldist(ax21_10,9.62)
plot_settings(ax21_10,time_Myr,t_ion,rho0str,alpha)
mcav = get_ionizedmass_inshell(x_src,y_src,z_src,0,0.53,'rad_21_10/new_nixyzhmf/nixyzhmf_00003.txt')
mchl_ion = get_ionizedmass_inshell(x_src,y_src,z_src,0.53,9.62,'rad_21_10/new_nixyzhmf/nixyzhmf_00003.txt')
mchl_neu = get_neutralmass_inshell(x_src,y_src,z_src,0.53,9.62,'rad_21_10/new_nixyzhmf/nixyzhmf_00003.txt')
print('21_10: ',mchl_ion,mchl_neu,mcav,mchl_neu/mcav)


##### rho 1e-20 #####

rho0str = r'$\rho_0 = 10^{-20}\ \mathrm{g\ cm^{-3}}$'
t_ff = 0.666 # Myr 

ax20_07 = axes[0,1]
alpha = 0.7
time0_Myr, sdf, sdf_sinks = read_dump('rad_20_07/cloud_20_07_clrsink1_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_20_07/cloud_20_07_clrsink1_21490') 
t_ion = time_Myr - time0_Myr 
plot_disttosrc_hist(ax20_07,x_src,y_src,z_src,'rad_20_07/new_nixyzhmf/nixyzhmf_00005.txt') 
plot_cavline(ax20_07,0.238)
plot_chnldist(ax20_07,6.54)
plot_settings(ax20_07,time_Myr,t_ion,rho0str,alpha)
mcav = get_ionizedmass_inshell(x_src,y_src,z_src,0,0.238,'rad_20_07/new_nixyzhmf/nixyzhmf_00005.txt')
mchl_ion = get_ionizedmass_inshell(x_src,y_src,z_src,0.238,6.54,'rad_20_07/new_nixyzhmf/nixyzhmf_00005.txt')
mchl_neu = get_neutralmass_inshell(x_src,y_src,z_src,0.238,6.54,'rad_20_07/new_nixyzhmf/nixyzhmf_00005.txt')
print('20_07: ',mchl_ion,mchl_neu,mcav,mchl_neu/mcav)


ax20_10 = axes[1,1]
alpha = 1.0
time0_Myr, sdf, sdf_sinks = read_dump('rad_20_10/cloud_20_10_clrsink14_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_20_10/cloud_20_10_clrsink14_20472') 
t_ion = time_Myr - time0_Myr 
plot_disttosrc_hist(ax20_10,x_src,y_src,z_src,'rad_20_10/new_nixyzhmf/nixyzhmf_00020.txt') 
plot_cavline(ax20_10,0.268)
plot_chnldist(ax20_10,5.96)
plot_settings(ax20_10,time_Myr,t_ion,rho0str,alpha)
mcav = get_ionizedmass_inshell(x_src,y_src,z_src,0,0.268,'rad_20_10/new_nixyzhmf/nixyzhmf_00020.txt')
mchl_ion = get_ionizedmass_inshell(x_src,y_src,z_src,0.268,5.96,'rad_20_10/new_nixyzhmf/nixyzhmf_00020.txt')
mchl_neu = get_neutralmass_inshell(x_src,y_src,z_src,0.268,5.96,'rad_20_10/new_nixyzhmf/nixyzhmf_00020.txt')
print('20_10: ',mchl_ion,mchl_neu,mcav,mchl_neu/mcav)



fig.tight_layout()
fig.subplots_adjust(hspace=0.03)

plt.savefig('chnl_ionization_intermed.png',dpi=200)



'''
Without channels 
'''

fig, axes = plt.subplots(nrows=2, ncols=2, subplot_kw=dict(frameon=True), figsize=(7,6))


##### rho 1e-21 #####

rho0str = r'$\rho_0 = 10^{-21}\ \mathrm{g\ cm^{-3}}$'
t_ff = 2.107 # Myr 

ax21_07 = axes[0,0]
alpha = 0.7
time0_Myr, sdf, sdf_sinks = read_dump('rad_21_07/cloud_21_07_clrsink4_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_21_07/cloud_21_07_clrsink4_29329')  
t_ion = time_Myr - time0_Myr 
plot_disttosrc_hist(ax21_07,x_src,y_src,z_src,'rad_21_07/new_nixyzhmf/nixyzhmf_00177.txt')
plot_cavline(ax21_07,2.97)
plot_chnldist(ax21_07,6.13)
plot_settings(ax21_07,time_Myr,t_ion,rho0str,alpha)
mcav = get_ionizedmass_inshell(x_src,y_src,z_src,0,2.97,'rad_21_07/new_nixyzhmf/nixyzhmf_00177.txt')
mchl_ion = get_ionizedmass_inshell(x_src,y_src,z_src,2.97,6.13,'rad_21_07/new_nixyzhmf/nixyzhmf_00177.txt')
mchl_neu = get_neutralmass_inshell(x_src,y_src,z_src,2.97,6.13,'rad_21_07/new_nixyzhmf/nixyzhmf_00177.txt')
print('21_07: ',mchl_ion,mchl_neu,mcav,mchl_neu/mcav)


ax21_10 = axes[1,0]
alpha = 1.0
time0_Myr, sdf, sdf_sinks = read_dump('rad_21_10/cloud_21_10_clrsink6_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_21_10/cloud_21_10_clrsink6_26245')  
t_ion = time_Myr - time0_Myr 
plot_disttosrc_hist(ax21_10,x_src,y_src,z_src,'rad_21_10/new_nixyzhmf/nixyzhmf_00098.txt') 
plot_cavline(ax21_10,4.18)
plot_chnldist(ax21_10,7.3)
plot_settings(ax21_10,time_Myr,t_ion,rho0str,alpha)
mcav = get_ionizedmass_inshell(x_src,y_src,z_src,0,4.18,'rad_21_10/new_nixyzhmf/nixyzhmf_00098.txt')
mchl_ion = get_ionizedmass_inshell(x_src,y_src,z_src,4.18,7.3,'rad_21_10/new_nixyzhmf/nixyzhmf_00098.txt')
mchl_neu = get_neutralmass_inshell(x_src,y_src,z_src,4.18,7.3,'rad_21_10/new_nixyzhmf/nixyzhmf_00098.txt')
print('21_10: ',mchl_ion,mchl_neu,mcav,mchl_neu/mcav)



##### rho 1e-20 #####

rho0str = r'$\rho_0 = 10^{-20}\ \mathrm{g\ cm^{-3}}$'
t_ff = 0.666 # Myr 

ax20_07 = axes[0,1]
alpha = 0.7
time0_Myr, sdf, sdf_sinks = read_dump('rad_20_07/cloud_20_07_clrsink1_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_20_07/cloud_20_07_clrsink1_24172') 
t_ion = time_Myr - time0_Myr 
plot_disttosrc_hist(ax20_07,x_src,y_src,z_src,'rad_20_07/new_nixyzhmf/nixyzhmf_00435.txt') 
plot_cavline(ax20_07,1.54)
plot_chnldist(ax20_07,4.62)
plot_settings(ax20_07,time_Myr,t_ion,rho0str,alpha)
mcav = get_ionizedmass_inshell(x_src,y_src,z_src,0,1.54,'rad_20_07/new_nixyzhmf/nixyzhmf_00435.txt')
mchl_ion = get_ionizedmass_inshell(x_src,y_src,z_src,1.54,4.62,'rad_20_07/new_nixyzhmf/nixyzhmf_00435.txt')
mchl_neu = get_neutralmass_inshell(x_src,y_src,z_src,1.54,4.62,'rad_20_07/new_nixyzhmf/nixyzhmf_00435.txt')
print('20_07: ',mchl_ion,mchl_neu,mcav,mchl_neu/mcav)


ax20_10 = axes[1,1]
alpha = 1.0
time0_Myr, sdf, sdf_sinks = read_dump('rad_20_10/cloud_20_10_clrsink14_00000')
isink, x_src, y_src, z_src = get_most_massive_sink(sdf_sinks)
time_Myr, sdf, sdf_sinks = read_dump('rad_20_10/cloud_20_10_clrsink14_22955') 
t_ion = time_Myr - time0_Myr 
plot_disttosrc_hist(ax20_10,x_src,y_src,z_src,'rad_20_10/new_nixyzhmf/nixyzhmf_00396.txt') 
plot_cavline(ax20_10,1.44)
plot_chnldist(ax20_10,4.7)
plot_settings(ax20_10,time_Myr,t_ion,rho0str,alpha)
mcav = get_ionizedmass_inshell(x_src,y_src,z_src,0,1.44,'rad_20_10/new_nixyzhmf/nixyzhmf_00396.txt')
mchl_ion = get_ionizedmass_inshell(x_src,y_src,z_src,1.44,4.7,'rad_20_10/new_nixyzhmf/nixyzhmf_00396.txt')
mchl_neu = get_neutralmass_inshell(x_src,y_src,z_src,1.44,4.7,'rad_20_10/new_nixyzhmf/nixyzhmf_00396.txt')
print('20_10: ',mchl_ion,mchl_neu,mcav,mchl_neu/mcav)




fig.tight_layout()
fig.subplots_adjust(hspace=0.03)

plt.savefig('chnl_ionization_final.png',dpi=200)


plt.show()



