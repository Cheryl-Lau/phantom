
import numpy as np
import matplotlib.pyplot as plt 
import sarracen
import glob
from scipy.spatial.transform import Rotation as R

import matplotlib
matplotlib.use("AGG")


utime = 4.706e14
udist = 3.086e+19
umass = 1.989e+33
unit_density = 6.768e-23
unit_velocity = 6.558e+3
unit_ergg = 4.301e+7
yr = 365*24*60*60 

t_ion = 1.2991e-1*utime/yr


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
    time_yr = time/yr
    sdf['rho'] = sdf['rho']*unit_density 
    sdf['u'] = sdf['u']*unit_ergg

    return time_yr, sdf, sdf_sinks



def plot_render_rho(ax,time_yr,sdf,xmin,xmax,ymin,ymax,vmin_in,vmax_in,viewangle):
    '''
    Plot input sdf log col-dens on the given ax
    '''
    ax = sdf.render('rho', ax=ax, xlim=(xmin,xmax), ylim=(ymin,ymax), rotation=[0,viewangle,0], rot_origin=[x_src,y_src,z_src], log_scale=True, cmap='gist_heat', vmin=vmin_in, vmax=vmax_in,  cbar_kws=dict(label=r'column $\rho$ [g/$\mathrm{cm}^2$]',orientation='vertical',shrink=0.7,pad=0.02))
    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.9*(ymax-ymin),r'$\mathrm{t_{ion}}$ = '+str(abs(round(time_yr,3)))+' yr', color='white')

    return 



def plot_render_u(ax,time_yr,sdf,xmin,xmax,ymin,ymax,vmin_in,vmax_in,viewangle):
    '''
    Plot input sdf log col-u on the given ax 
    '''
    ax = sdf.render('u', ax=ax, xlim=(xmin,xmax), ylim=(ymin,ymax), rotation=[0,viewangle,0], rot_origin=[x_src,y_src,z_src], log_scale=True, cmap='gnuplot2', vmin=vmin_in, vmax=vmax_in,  cbar_kws=dict(label=r'column $\rho$-weighted $u$ [$\mathrm{erg/g\ cm}$]',orientation='vertical',shrink=0.7,pad=0.02), dens_weight=True)
    ax.set_xlabel('x [pc]')
    ax.set_ylabel('y [pc]')
    ax.text(xmin + 0.1*(xmax-xmin), ymin + 0.9*(ymax-ymin), r'$\mathrm{t_{ion}}$ = '+str(abs(round(time_yr,3)))+' yr', color='white')

    return 



def plot_sinks(ax,sdf_sinks,viewangle): 

    xyz_sink = np.dstack([sdf_sinks['x'].values,sdf_sinks['y'].values,sdf_sinks['z'].values])[0]

    r = R.from_euler('zyx', [0,viewangle,0], degrees=True)

    xyz_proj = r.apply(xyz_sink)
    xyz_proj = xyz_proj.transpose() 

    ax.scatter(xyz_proj[0],xyz_proj[1],s=0.1,color='white')

    return 




rho0str = r'$\rho_0 = 10^{-21}\ \mathrm{g\ cm^{-3}}$'

x_src = 0
y_src = 0
z_src = 0

radius = 6.5
xmin = x_src - radius 
xmax = x_src + radius 
ymin = y_src - radius 
ymax = y_src + radius 

rho_min = 1e+0
rho_max = 1e+4

u_min = 1e+10
u_max = 1e+14

plt.gca().set_aspect('equal')

viewangle = 0.1



for dumpfile in sorted(glob.glob('./cloud_21_07_clrsink4_*00')):  
    print(dumpfile)
    time_yr, sdf, sdf_sinks = read_dump(dumpfile)

    viewangle += 0.3


    # Column Density

    fig1, ax1 = plt.subplots(nrows=1, ncols=1, subplot_kw=dict(frameon=True), figsize=(8,8))
    plot_render_rho(ax1,time_yr-t_ion,sdf,xmin,xmax,ymin,ymax,rho_min,rho_max,viewangle)
    plot_sinks(ax1,sdf_sinks,viewangle)

    imfilename = 'rhorender_'+dumpfile[-5:]+'.png'
    fig1.savefig(imfilename,dpi=200)


    # Column Internal energy 

    fig2, ax2 = plt.subplots(nrows=1, ncols=1, subplot_kw=dict(frameon=True), figsize=(8,8))
    plot_render_u(ax2,time_yr-t_ion,sdf,xmin,xmax,ymin,ymax,u_min,u_max,viewangle)
    plot_sinks(ax2,sdf_sinks,viewangle)

    imfilename = 'urender_'+dumpfile[-5:]+'.png'
    fig2.savefig(imfilename,dpi=200)











