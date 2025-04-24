
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib._layoutgrid import plot_children
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



def plot_render_u(ax,time_Myr,sdf,xmin,xmax,ymin,ymax,vmin_in,vmax_in):
    '''
    Plot input sdf log col-u on the given ax 
    '''
    ax = sdf.render('u', ax=ax, xlim=(xmin,xmax), ylim=(ymin,ymax), vmin=vmin_in,vmax=vmax_in, log_scale=True, cmap='gnuplot2', cbar_kws=dict(label=r'column $\rho$-weighted $u$ [$\mathrm{erg/g\ cm}$]',orientation='vertical',shrink=0.7,pad=0.01), dens_weight=True)
    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.9*(ymax-ymin), r'$t_\mathrm{SN} = $'+str(round(time_Myr,2))+' Myr', color='white')

    return 



fig = plt.figure(figsize=(22,7),layout="constrained")
gs = fig.add_gridspec(2, 5)
ax_semiconf = fig.add_subplot(gs[:, :2])
ax_ffcloud = fig.add_subplot(gs[0, 2])
ax_ffenv = fig.add_subplot(gs[0, 3])
ax_ffsmear = fig.add_subplot(gs[0, 4])
ax_tffcloud = fig.add_subplot(gs[1, 2])
ax_tffenv = fig.add_subplot(gs[1, 3])
ax_tffsmear = fig.add_subplot(gs[1, 4])


xmin = -10
xmax = 10
ymin = -10
ymax = 10


# Semi-confined SN 
u_min = 1e+9
u_max = 1e+13

time_SN = 0.486
time_Myr, sdf = read_dump('turbcloud_semiconf/cloud_20_10_clrsink14_04500')
plot_render_u(ax_semiconf,time_Myr-time_SN,sdf,xmin,xmax,ymin,ymax,u_min,u_max)
ax_semiconf.set_title('Semi-confined SN')


# ff_cloud 
u_min = 1e+12
u_max = 7e+12

time_SN = 50.002
time_Myr, sdf = read_dump('ff_cloud/ffcloud_02000')
plot_render_u(ax_ffcloud,time_Myr-time_SN,sdf,xmin,xmax,ymin,ymax,u_min,u_max)
ax_ffcloud.set_title(r'FF - $\rho_\mathrm{cloud}$')

time_SN = 1.200
time_Myr, sdf = read_dump('turbff_cloud/turbffcloud_01900')
plot_render_u(ax_tffcloud,time_Myr-time_SN,sdf,xmin,xmax,ymin,ymax,u_min,u_max)
ax_tffcloud.set_title(r'Turb FF - $\rho_\mathrm{cloud}$')


# ff_env 
u_min = 6e+13
u_max = 4e+14

time_SN = 6.000
time_Myr, sdf = read_dump('ff_env/ffenv_01900')
plot_render_u(ax_ffenv,time_Myr-time_SN,sdf,xmin,xmax,ymin,ymax,u_min,u_max)
ax_ffenv.set_title(r'FF - $\rho_\mathrm{env}$')

time_SN = 1.200
time_Myr, sdf = read_dump('turbff_env/turbffenv_01900')
plot_render_u(ax_tffenv,time_Myr-time_SN,sdf,xmin,xmax,ymin,ymax,u_min,u_max)
ax_tffenv.set_title(r'Turb FF - $\rho_\mathrm{env}$')


# ff_smear 
u_min = 1e+10
u_max = 1e+16

time_SN = 60.002
time_Myr, sdf = read_dump('ff_smear/ffsmear_02000')
plot_render_u(ax_ffsmear,time_Myr-time_SN,sdf,xmin,xmax,ymin,ymax,u_min,u_max)
ax_ffsmear.set_title(r'FF - $\rho_\mathrm{smear}$')

time_SN = 1.500
time_Myr, sdf = read_dump('turbff_smear/turbffsmear_02000')
plot_render_u(ax_tffsmear,time_Myr-time_SN,sdf,xmin,xmax,ymin,ymax,u_min,u_max)
ax_tffsmear.set_title(r'Turb FF - $\rho_\mathrm{smear}$')


#plot_children(fig)

plt.savefig('turbboxes.png',dpi=200)
plt.savefig('turbboxes.pdf',format='pdf')
plt.show()























