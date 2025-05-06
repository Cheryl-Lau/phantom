'''
Plot snapshots of energy and momentum across longitude 
'''

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
solarm = 1.989e+33
pc = 3.086E+18
km = 1e5
mass_proton_cgs = 1.67262158e-24


# Time range 
tmin = 0.0
tmax = 0.059


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



def plot_line(ax,x,y,colour,labelstr):

    x_sorted = sorted(x)
    y_sorted = [ydum for _, ydum in sorted(zip(x,y))]
    x_sorted.pop(0)
    y_sorted.pop(0)
    ax.set_yscale('log')
    ax.plot(x_sorted,y_sorted,color=colour,linewidth=1,label=labelstr)    
    ax.grid('on', linestyle='--')

    ax.xaxis.set_tick_params(which='both', labelbottom=True)
    ax.yaxis.set_tick_params(which='both', labelbottom=False)

    return 


def plot_dashed_line(ax,x,y,colour,labelstr):

    x_sorted = sorted(x)
    y_sorted = [ydum for _, ydum in sorted(zip(x,y))]
    x_sorted.pop(0)
    y_sorted.pop(0)
    ax.set_yscale('log')
    ax.plot(x_sorted,y_sorted,'--',color=colour,linewidth=1,label=labelstr)    

    return 


def plot_dotted_line(ax,x,y,colour,labelstr):

    x_sorted = sorted(x)
    y_sorted = [ydum for _, ydum in sorted(zip(x,y))]
    x_sorted.pop(0)
    y_sorted.pop(0)
    ax.set_yscale('log')
    ax.plot(x_sorted,y_sorted,':',color=colour,linewidth=1,label=labelstr)    

    return 



def read_energy_momentum(filepath,dumpfile):


    KE = np.loadtxt(filepath+'/theta_energies/theta_ekin_'+dumpfile+'.dat')
    TE = np.loadtxt(filepath+'/theta_energies/theta_etherm_'+dumpfile+'.dat')
    TotE = KE + TE 
    p = np.loadtxt(filepath+'/theta_energies/theta_momen_'+dumpfile+'.dat')

    ntheta = len(KE)
    theta = np.linspace(-180,180,ntheta)

    return theta, TotE, KE, TE, p 



fig, axes = plt.subplots(nrows=4, ncols=4, subplot_kw=dict(frameon=True), figsize=(10,8))
ax1_TotE = axes[0,0]
ax2_TotE = axes[0,1]
ax3_TotE = axes[0,2]
ax4_TotE = axes[0,3]
ax1_KE = axes[1,0]
ax2_KE = axes[1,1]
ax3_KE = axes[1,2]
ax4_KE = axes[1,3]
ax1_TE = axes[2,0]
ax2_TE = axes[2,1]
ax3_TE = axes[2,2]
ax4_TE = axes[2,3]
ax1_p = axes[3,0]
ax2_p = axes[3,1]
ax3_p = axes[3,2]
ax4_p = axes[3,3]


xmin = -180
xmax = 180
ymin_energ = 1e+45
ymax_energ = 3e+49
ymin_momen = 1e+39
ymax_momen = 3.5e+41

axes_TotE = [ax1_TotE,ax2_TotE,ax3_TotE,ax4_TotE]
for ax in axes_TotE:
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin_energ,ymax_energ])
    ax.axhline(y=1e51/360, color='black', linestyle='--',label=r'$E_\mathrm{tot} = 10^{51}/360$')
ax1_TotE.set_ylabel(r'$E_\mathrm{tot}$ [erg]')

axes_KE = [ax1_KE,ax2_KE,ax3_KE,ax4_KE]
for ax in axes_KE:
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin_energ,ymax_energ])
ax1_KE.set_ylabel(r'$E_\mathrm{kin}$ [erg]')

ax_TE = [ax1_TE,ax2_TE,ax3_TE,ax4_TE]
for ax in ax_TE:
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin_energ,ymax_energ])
ax1_TE.set_ylabel(r'$E_\mathrm{therm}$ [erg]')

axes_momen = [ax1_p,ax2_p,ax3_p,ax4_p]
for ax in axes_momen:
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin_momen,ymax_momen])
    ax.set_xlabel(r'$\theta$ [$^{\circ}$]')
ax1_p.set_ylabel(r'$\mu_r$ [$\mathrm{g\ cm\ s^{-1}}$]')


axes_t1 = [ax1_TotE,ax1_KE,ax1_TE]
for ax in axes_t1:
    ax.text(-170,1e+49,'t = 0.0005 Myr')
ax1_p.text(-170,2e+41,'t = 0.0005 Myr')

axes_t2 = [ax2_TotE,ax2_KE,ax2_TE]
for ax in axes_t2:
    ax.text(-170,1e+49,'t = 0.001 Myr')
ax2_p.text(-170,2e+41,'t = 0.001 Myr')

axes_t3 = [ax3_TotE,ax3_KE,ax3_TE]
for ax in axes_t3:
    ax.text(-170,1e+49,'t = 0.01 Myr')
ax3_p.text(-170,2e+41,'t = 0.01 Myr')

axes_t4 = [ax4_TotE,ax4_KE,ax4_TE]
for ax in axes_t4:
    ax.text(-170,1e+49,'t = 0.06 Myr')
ax4_p.text(-170,2e+41,'t = 0.06 Myr')



#
# FF env 
#

filepath = 'ff_env4'

dumpfile = 'freefield_00000'
time_Myr0, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE0, KE0, TE0, p0 = read_energy_momentum(filepath,dumpfile)

dumpfile = 'freefield_00000'  #'freefield_00010'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_line(ax1_TotE,theta,TotE-TotE0,'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
plot_line(ax1_KE,theta,KE-KE0,'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
plot_line(ax1_TE,theta,TE-TE0,'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
plot_line(ax1_p,theta,abs(p-p0),'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
print('ffenv t1',time_Myr-time_Myr0)

dumpfile = 'freefield_00019'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_line(ax2_TotE,theta,TotE-TotE0,'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
plot_line(ax2_KE,theta,KE-KE0,'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
plot_line(ax2_TE,theta,TE-TE0,'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
plot_line(ax2_p,theta,abs(p-p0),'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
print('ffenv t2',time_Myr-time_Myr0)

dumpfile = 'freefield_00190'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_line(ax3_TotE,theta,TotE-TotE0,'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
plot_line(ax3_KE,theta,KE-KE0,'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
plot_line(ax3_TE,theta,TE-TE0,'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
plot_line(ax3_p,theta,abs(p-p0),'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
print('ffenv t3',time_Myr-time_Myr0)

dumpfile = 'freefield_01150'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_line(ax4_TotE,theta,TotE-TotE0,'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
plot_line(ax4_KE,theta,KE-KE0,'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
plot_line(ax4_TE,theta,TE-TE0,'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
plot_line(ax4_p,theta,abs(p-p0),'cornflowerblue',r'FF - $\rho_\mathrm{env}$')
print('ffenv t4',time_Myr-time_Myr0)


#
# Turb FF env 
#

filepath = 'turbff_env4'

dumpfile = 'turbffenv_00000'
time_Myr0, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE0, KE0, TE0, p0 = read_energy_momentum(filepath,dumpfile)

dumpfile = 'turbffenv_00010'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_dotted_line(ax1_TotE,theta,TotE-TotE0,'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
plot_dotted_line(ax1_KE,theta,KE-KE0,'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
plot_dotted_line(ax1_TE,theta,TE-TE0,'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
plot_dotted_line(ax1_p,theta,abs(p-p0),'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
print('tffenv t1',time_Myr-time_Myr0)

dumpfile = 'turbffenv_00019'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_dotted_line(ax2_TotE,theta,TotE-TotE0,'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
plot_dotted_line(ax2_KE,theta,KE-KE0,'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
plot_dotted_line(ax2_TE,theta,TE-TE0,'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
plot_dotted_line(ax2_p,theta,abs(p-p0),'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
print('tffenv t2',time_Myr-time_Myr0)

dumpfile = 'turbffenv_00190'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_dotted_line(ax3_TotE,theta,TotE-TotE0,'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
plot_dotted_line(ax3_KE,theta,KE-KE0,'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
plot_dotted_line(ax3_TE,theta,TE-TE0,'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
plot_dotted_line(ax3_p,theta,abs(p-p0),'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
print('tffenv t3',time_Myr-time_Myr0)

dumpfile = 'turbffenv_01150'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_dotted_line(ax4_TotE,theta,TotE-TotE0,'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
plot_dotted_line(ax4_KE,theta,KE-KE0,'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
plot_dotted_line(ax4_TE,theta,TE-TE0,'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
plot_dotted_line(ax4_p,theta,abs(p-p0),'cornflowerblue',r'TurbFF - $\rho_\mathrm{env}$')
print('tffenv t4',time_Myr-time_Myr0)


#
# FF smear 
#

filepath = 'ff_smear2'

dumpfile = 'ffsmear_00000'
time_Myr0, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE0, KE0, TE0, p0 = read_energy_momentum(filepath,dumpfile)

dumpfile = 'ffsmear_00010'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_line(ax1_TotE,theta,TotE-TotE0,'navy',r'FF - $\rho_\mathrm{smear}$')
plot_line(ax1_KE,theta,KE-KE0,'navy',r'FF - $\rho_\mathrm{smear}$')
plot_line(ax1_TE,theta,TE-TE0,'navy',r'FF - $\rho_\mathrm{smear}$')
plot_line(ax1_p,theta,abs(p-p0),'navy',r'FF - $\rho_\mathrm{smear}$')
print('ffsmear t1',time_Myr-time_Myr0)

dumpfile = 'ffsmear_00019'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_line(ax2_TotE,theta,TotE-TotE0,'navy',r'FF - $\rho_\mathrm{smear}$')
plot_line(ax2_KE,theta,KE-KE0,'navy',r'FF - $\rho_\mathrm{smear}$')
plot_line(ax2_TE,theta,TE-TE0,'navy',r'FF - $\rho_\mathrm{smear}$')
plot_line(ax2_p,theta,abs(p-p0),'navy',r'FF - $\rho_\mathrm{smear}$')
print('ffsmear t2',time_Myr-time_Myr0)

dumpfile = 'ffsmear_00190'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_line(ax3_TotE,theta,TotE-TotE0,'navy',r'FF - $\rho_\mathrm{smear}$')
plot_line(ax3_KE,theta,KE-KE0,'navy',r'FF - $\rho_\mathrm{smear}$')
plot_line(ax3_TE,theta,TE-TE0,'navy',r'FF - $\rho_\mathrm{smear}$')
plot_line(ax3_p,theta,abs(p-p0),'navy',r'FF - $\rho_\mathrm{smear}$')
print('ffsmear t3',time_Myr-time_Myr0)

dumpfile = 'ffsmear_01150'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_line(ax4_TotE,theta,TotE-TotE0,'navy',r'FF - $\rho_\mathrm{smear}$')
plot_line(ax4_KE,theta,KE-KE0,'navy',r'FF - $\rho_\mathrm{smear}$')
plot_line(ax4_TE,theta,TE-TE0,'navy',r'FF - $\rho_\mathrm{smear}$')
plot_line(ax4_p,theta,abs(p-p0),'navy',r'FF - $\rho_\mathrm{smear}$')
print('ffsmear t4',time_Myr-time_Myr0)


#
# Turb FF smear 
#

filepath = 'turbff_smear'

dumpfile = 'turbffsmear_00000'
time_Myr0, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE0, KE0, TE0, p0 = read_energy_momentum(filepath,dumpfile)

dumpfile = 'turbffsmear_00100'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_dotted_line(ax1_TotE,theta,TotE-TotE0,'navy',r'TurbFF - $\rho_\mathrm{smear}$')
plot_dotted_line(ax1_KE,theta,KE-KE0,'navy',r'TurbFF - $\rho_\mathrm{smear}$')
plot_dotted_line(ax1_TE,theta,TE-TE0,'navy',r'TurbFF - $\rho_\mathrm{smear}$')
plot_dotted_line(ax1_p,theta,abs(p-p0),'navy',r'TurbFF - $\rho_\mathrm{smear}$')
print('tffsmear t1',time_Myr-time_Myr0)

dumpfile = 'turbffsmear_00190'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_dotted_line(ax2_TotE,theta,TotE-TotE0,'navy',r'TurbFF - $\rho_\mathrm{smear}$')
plot_dotted_line(ax2_KE,theta,KE-KE0,'navy',r'TurbFF - $\rho_\mathrm{smear}$')
plot_dotted_line(ax2_TE,theta,TE-TE0,'navy',r'TurbFF - $\rho_\mathrm{smear}$')
plot_dotted_line(ax2_p,theta,abs(p-p0),'navy',r'TurbFF - $\rho_\mathrm{smear}$')
print('tffsmear t2',time_Myr-time_Myr0)

dumpfile = 'turbffsmear_01900'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_dotted_line(ax3_TotE,theta,TotE-TotE0,'navy',r'TurbFF - $\rho_\mathrm{smear}$')
plot_dotted_line(ax3_KE,theta,KE-KE0,'navy',r'TurbFF - $\rho_\mathrm{smear}$')
plot_dotted_line(ax3_TE,theta,TE-TE0,'navy',r'TurbFF - $\rho_\mathrm{smear}$')
plot_dotted_line(ax3_p,theta,abs(p-p0),'navy',r'TurbFF - $\rho_\mathrm{smear}$')
print('tffsmear t3',time_Myr-time_Myr0)

dumpfile = 'turbffsmear_11500'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_dotted_line(ax4_TotE,theta,TotE-TotE0,'navy',r'TurbFF - $\rho_\mathrm{smear}$')
plot_dotted_line(ax4_KE,theta,KE-KE0,'navy',r'TurbFF - $\rho_\mathrm{smear}$')
plot_dotted_line(ax4_TE,theta,TE-TE0,'navy',r'TurbFF - $\rho_\mathrm{smear}$')
plot_dotted_line(ax4_p,theta,abs(p-p0),'navy',r'TurbFF - $\rho_\mathrm{smear}$')
print('tffsmear t4',time_Myr-time_Myr0)



#
# Cloud 
#

filepath = 'turbcloud_semiconf'

dumpfile = 'cloud_20_10_clrsink14_01350'
time_Myr0, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE0, KE0, TE0, p0 = read_energy_momentum(filepath,dumpfile)

dumpfile = 'cloud_20_10_clrsink14_03455'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_line(ax1_TotE,theta,TotE-TotE0,'red','cloud')
plot_line(ax1_KE,theta,KE-KE0,'red','cloud')
plot_line(ax1_TE,theta,TE-TE0,'red','cloud')
plot_line(ax1_p,theta,abs(p-p0),'red','cloud')
print('cloud t1',time_Myr-time_Myr0)

dumpfile = 'cloud_20_10_clrsink14_03950'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_line(ax2_TotE,theta,TotE-TotE0,'red','cloud')
plot_line(ax2_KE,theta,KE-KE0,'red','cloud')
plot_line(ax2_TE,theta,TE-TE0,'red','cloud')
plot_line(ax2_p,theta,abs(p-p0),'red','cloud')
print('cloud t2',time_Myr-time_Myr0)

dumpfile = 'cloud_20_10_clrsink14_04550'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_line(ax3_TotE,theta,TotE-TotE0,'red','cloud')
plot_line(ax3_KE,theta,KE-KE0,'red','cloud')
plot_line(ax3_TE,theta,TE-TE0,'red','cloud')
plot_line(ax3_p,theta,abs(p-p0),'red','cloud')
print('cloud t3',time_Myr-time_Myr0)

dumpfile = 'cloud_20_10_clrsink14_09730'
time_Myr, sdf = read_dump(filepath+'/'+dumpfile)
theta, TotE, KE, TE, p = read_energy_momentum(filepath,dumpfile)
plot_line(ax4_TotE,theta,TotE-TotE0,'red','cloud')
plot_line(ax4_KE,theta,KE-KE0,'red','cloud')
plot_line(ax4_TE,theta,TE-TE0,'red','cloud')
plot_line(ax4_p,theta,abs(p-p0),'red','cloud')
print('cloud t4',time_Myr-time_Myr0)



for ax in axes_t1:
    ax.yaxis.set_tick_params(which='both', labelbottom=True)
ax1_p.yaxis.set_tick_params(which='both', labelbottom=True)

for ax in axes_t4:
    ax.legend(loc='lower right',fontsize=7)
ax4_p.legend(loc='lower right',fontsize=7)


fig.tight_layout(pad=0.2)


fig.savefig('energy_momentum_theta.pdf',format='pdf')
fig.savefig('energy_momentum_theta.png')

plt.show()


















