
'''
Plot the evolution of energy and momentum contained within spheres of different radii
'''

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import glob 
import os 



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



def plot_line(ax,x,y,colour,labelstr):

    x_sorted = sorted(x)
    y_sorted = [ydum for _, ydum in sorted(zip(x,y))]
    x_sorted.pop(0)
    y_sorted.pop(0)
    ax.plot(x_sorted,y_sorted,color=colour,linewidth=1,label=labelstr)    
    ax.set_yscale('log')
    ax.set_xlim([tmin,tmax])
    ax.grid('on', linestyle='--')

    ax.xaxis.set_tick_params(which='both', labelbottom=True)
    ax.yaxis.set_tick_params(which='both', labelbottom=False)

    return 


def plot_dashed_line(ax,x,y,colour,labelstr):

    x_sorted = sorted(x)
    y_sorted = [ydum for _, ydum in sorted(zip(x,y))]
    x_sorted.pop(0)
    y_sorted.pop(0)
    ax.plot(x_sorted,y_sorted,'--',color=colour,linewidth=1,label=labelstr)    
    ax.set_yscale('log')
    ax.set_xlim([tmin,tmax])

    return 


def plot_dotted_line(ax,x,y,colour,labelstr):

    x_sorted = sorted(x)
    y_sorted = [ydum for _, ydum in sorted(zip(x,y))]
    x_sorted.pop(0)
    y_sorted.pop(0)
    ax.plot(x_sorted,y_sorted,':',color=colour,linewidth=2,label=labelstr)    
    ax.set_yscale('log')
    ax.set_xlim([tmin,tmax])

    return 


def get_cioffi(rho_cgs):

    nH = rho_cgs/mass_proton_cgs
    momen_cioffi88 = 4.8e5*km * (1e0)**(13e0/14e0) * (1e0)**(-3e0/14e0) * (nH)**(-1e0/7e0) * solarm

    return momen_cioffi88


def get_kimmcen(rho_cgs):

    nH = rho_cgs/mass_proton_cgs
    momen_kimmcen14 = 3e5*km * (1e0)**(16e0/17e0) * (nH)**(-2e0/17e0) * solarm * (1e0)**(-0.14)

    return momen_kimmcen14



def read_rawfiles(filepath_energ,filepath_momen,sim_label,rho0_ff,df_energ,df_momen):

    # Energy files 
    first = True 
    for filename in sorted(glob.glob(filepath_energ)):  # over_time 
        time = np.loadtxt(filename,skiprows=1,max_rows=1,unpack=True)
        radius, ke, te, tote = np.loadtxt(filename,skiprows=3,unpack=True)
        if (first == True):
            time0 = time
            ke0 = ke
            te0 = te
            tote0 = tote
            first = False 
        df_add = pd.DataFrame({'sim':sim_label,'time':(time-time0)/Myr,'radius':radius/pc,'ekin':ke-ke0,'etherm':te-te0,'etot':tote-tote0})
        df_energ = pd.concat([df_energ,df_add],ignore_index=True)

    # Momentum files 
    first = True
    for filename in sorted(glob.glob(filepath_momen)):  # over_time 
        with open(filename, 'r') as file:
            lines = file.readlines()
            time = float(lines[1].strip())
            mom_cioffi = get_cioffi(rho0_ff)
            mom_kimmcen = get_kimmcen(rho0_ff)
        radius, mom = np.loadtxt(filename,skiprows=6,unpack=True)
        if (first == True):
            time0 = time
            mom0 = mom
            first = False 
        df_add = pd.DataFrame({'sim':sim_label,'time':(time-time0)/Myr,'radius':radius/pc,'momen':abs(mom-mom0),'momen_cioffi':mom_cioffi,'momen_kimmcen':mom_kimmcen})
        df_momen = pd.concat([df_momen,df_add],ignore_index=True)

    # Store spheres radii info 
    global nrad, radii
    nrad = len(radius)
    radii = radius/pc

    return df_energ, df_momen



def combine_frames():

    df_energ = pd.DataFrame(columns=('sim', 'time', 'radius', 'ekin', 'etherm', 'etot'))
    df_momen = pd.DataFrame(columns=('sim', 'time', 'radius', 'momen', 'momen_cioffi', 'momen_kimmcen'))

    # Turb cloud 
    filepath_energ = 'turbcloud_semiconf/energy_insphere/energy_insphere_cloud_20_10_clrsink14_*.dat'
    filepath_momen = 'turbcloud_semiconf/momentum_insphere/momentum_insphere_cloud_20_10_clrsink14_*.dat'
    df_energ, df_momen = read_rawfiles(filepath_energ, filepath_momen, 'cloud', 1e-20, df_energ, df_momen)

    # FF_env 
    filepath_energ = 'ff_env4/energy_insphere/energy_insphere_freefield_*.dat'
    filepath_momen = 'ff_env4/momentum_insphere/momentum_insphere_freefield_*.dat'
    df_energ, df_momen = read_rawfiles(filepath_energ, filepath_momen, 'ffenv', 4e-25, df_energ, df_momen)

    # TurbFF_env 
    filepath_energ = 'turbff_env4/energy_insphere/energy_insphere_turbffenv_*.dat'
    filepath_momen = 'turbff_env4/momentum_insphere/momentum_insphere_turbffenv_*.dat'
    df_energ, df_momen = read_rawfiles(filepath_energ, filepath_momen, 'tffenv', 4e-25, df_energ, df_momen)

    # FF_smear 
    filepath_energ = 'ff_smear2/energy_insphere/energy_insphere_ffsmear_*.dat'
    filepath_momen = 'ff_smear2/momentum_insphere/momentum_insphere_ffsmear_*.dat'
    df_energ, df_momen = read_rawfiles(filepath_energ, filepath_momen, 'ffsmear', 8.46e-23, df_energ, df_momen)

    # TurbFF_smear 
    filepath_energ = 'turbff_smear/energy_insphere/energy_insphere_turbffsmear_*.dat'
    filepath_momen = 'turbff_smear/momentum_insphere/momentum_insphere_turbffsmear_*.dat'
    df_energ, df_momen = read_rawfiles(filepath_energ, filepath_momen, 'tffsmear', 8.46e-23, df_energ, df_momen)

    print(df_energ)
    print(df_momen)

    return df_energ, df_momen



def plot_energ(ax,rad_pc,energ_type,df_energ): 

    # Turb cloud 
    iselect = np.where(np.logical_and(df_energ['radius'] == rad_pc, df_energ['sim']=='cloud'))[0]
    time = df_energ['time'].iloc[iselect].to_numpy()
    energ = df_energ[energ_type].iloc[iselect].to_numpy()
    plot_line(ax,time,energ,'red','Cloud')

    # FF_env 
    iselect = np.where(np.logical_and(df_energ['radius'] == rad_pc, df_energ['sim']=='ffenv'))[0]
    time = df_energ['time'].iloc[iselect].to_numpy()
    energ = df_energ[energ_type].iloc[iselect].to_numpy()
    plot_line(ax,time,energ,'lightblue','FF env')

    # TurbFF_env 
    iselect = np.where(np.logical_and(df_energ['radius'] == rad_pc, df_energ['sim']=='tffenv'))[0]
    time = df_energ['time'].iloc[iselect].to_numpy()
    energ = df_energ[energ_type].iloc[iselect].to_numpy()
    plot_line(ax,time,energ,'royalblue','TurbFF env')

    # FF_smear 
    iselect = np.where(np.logical_and(df_energ['radius'] == rad_pc, df_energ['sim']=='ffsmear'))[0]
    time = df_energ['time'].iloc[iselect].to_numpy()
    energ = df_energ[energ_type].iloc[iselect].to_numpy()
    plot_dashed_line(ax,time,energ,'lightblue','FF smear')

    # TurbFF_smear 
    iselect = np.where(np.logical_and(df_energ['radius'] == rad_pc, df_energ['sim']=='tffsmear'))[0]
    time = df_energ['time'].iloc[iselect].to_numpy()
    energ = df_energ[energ_type].iloc[iselect].to_numpy()
    plot_dashed_line(ax,time,energ,'royalblue','TurbFF smear')

    return 


def plot_momen(ax,rad_pc,df_momen): 

    # Cloud 
    iselect = np.where(np.logical_and(df_momen['radius'] == rad_pc, df_momen['sim']=='cloud'))[0]
    time = df_momen['time'].iloc[iselect].to_numpy()
    momen = df_momen['momen'].iloc[iselect].to_numpy()
    plot_line(ax,time,momen,'red','Cloud')

    # FF_env 
    iselect = np.where(np.logical_and(df_momen['radius'] == rad_pc, df_momen['sim']=='ffenv'))[0]
    time = df_momen['time'].iloc[iselect].to_numpy()
    momen = df_momen['momen'].iloc[iselect].to_numpy()
    plot_line(ax,time,momen,'lightblue','FF env')

    momen_cioffi = df_momen['momen_cioffi'].iloc[iselect].to_numpy()
    momen_kimmcen = df_momen['momen_kimmcen'].iloc[iselect].to_numpy()
    plot_dotted_line(ax,time,momen_cioffi,'black',r'Cioffi88 for $\rho_\mathrm{env}$')
#    plot_dotted_line(ax,time,momen_kimmcen,'grey',r'KimmCen14 for $\rho_\mathrm{env}$')

    # TurbFF_env 
    iselect = np.where(np.logical_and(df_momen['radius'] == rad_pc, df_momen['sim']=='tffenv'))[0]
    time = df_momen['time'].iloc[iselect].to_numpy()
    momen = df_momen['momen'].iloc[iselect].to_numpy()
    plot_line(ax,time,momen,'royalblue','TurbFF env')

    # FF_smear 
    iselect = np.where(np.logical_and(df_momen['radius'] == rad_pc, df_momen['sim']=='ffsmear'))[0]
    time = df_momen['time'].iloc[iselect].to_numpy()
    momen = df_momen['momen'].iloc[iselect].to_numpy()
    plot_dashed_line(ax,time,momen,'lightblue','FF smear')

    momen_cioffi = df_momen['momen_cioffi'].iloc[iselect].to_numpy()
    momen_kimmcen = df_momen['momen_kimmcen'].iloc[iselect].to_numpy()
    plot_dotted_line(ax,time,momen_cioffi,'grey',r'Cioffi88 for $\rho_\mathrm{smear}$')
#    plot_dotted_line(ax,time,momen_kimmcen,'grey',r'KimmCen14 for $\rho_\mathrm{smear}$')

    # TurbFF_smear 
    iselect = np.where(np.logical_and(df_momen['radius'] == rad_pc, df_momen['sim']=='tffsmear'))[0]
    time = df_momen['time'].iloc[iselect].to_numpy()
    momen = df_momen['momen'].iloc[iselect].to_numpy()
    plot_dashed_line(ax,time,momen,'royalblue','TurbFF smear')

    return 


def plot_energ_for_all_rad(axes,energ_type,df_energ):

    if (len(axes) != nrad):
        print('number of radii and number of energy axes are not matching!')

    for irad,ax in enumerate(axes): 
        rad = radii[irad]
        plot_energ(ax,rad,energ_type,df_energ)

    return 


def plot_momen_for_all_rad(axes,df_momen):

    if (len(axes) != nrad):
        print('number of radii and number of momentum axes are not matching!')

    for irad,ax in enumerate(axes): 
        rad = radii[irad]
        plot_momen(ax,rad,df_momen)
        print(ax)

    return 



def main():

    df_energ, df_momen = combine_frames()

    fig = plt.figure(figsize=[9,10],dpi=100)

    # Total energy
    ax1 = fig.add_subplot(431)
    ax2 = fig.add_subplot(432)
    ax3 = fig.add_subplot(433)
    axes = [ax1,ax2,ax3]
    plot_energ_for_all_rad(axes,'etot',df_energ)
    for ax in axes:
        ax.set_ylim([1E46,1E52])
    ax1.set_ylabel('total energy [erg]')
    ax1.yaxis.set_tick_params(which='both', labelbottom=True)
    ax1.set_title('r = 4 pc')
    ax2.set_title('r = 7 pc')
    ax3.set_title('r = 10 pc')
    ax3.legend(fontsize=6.5,loc='lower right')

    # Kinetic energy 
    ax4 = fig.add_subplot(434)
    ax5 = fig.add_subplot(435)
    ax6 = fig.add_subplot(436)
    axes = [ax4,ax5,ax6]
    plot_energ_for_all_rad(axes,'ekin',df_energ)
    for ax in axes:
        ax.set_ylim([1E43,1E53])
    ax4.set_ylabel('kinetic energy [erg]')
    ax4.yaxis.set_tick_params(which='both', labelbottom=True)
    ax6.legend(fontsize=6.5,loc='lower right')

    # Thermal energy 
    ax7 = fig.add_subplot(437)
    ax8 = fig.add_subplot(438)
    ax9 = fig.add_subplot(439)
    axes = [ax7,ax8,ax9]
    plot_energ_for_all_rad(axes,'etherm',df_energ)
    for ax in axes:
        ax.set_ylim([1E45,1E52])
    ax7.set_ylabel('thermal energy [erg]')
    ax7.yaxis.set_tick_params(which='both', labelbottom=True)
    ax9.legend(fontsize=6.5,loc='lower right')

    # Momentum 
    ax10 = fig.add_subplot(4,3,10)
    ax11 = fig.add_subplot(4,3,11)
    ax12 = fig.add_subplot(4,3,12)
    axes = [ax10,ax11,ax12]
    plot_momen_for_all_rad(axes,df_momen)
    for ax in axes:
        ax.set_ylim([1E36,1E46])
    ax10.set_xlabel('time [Myr]')
    ax11.set_xlabel('time [Myr]')
    ax12.set_xlabel('time [Myr]')
    ax10.set_ylabel('radial momentum [$\mathrm{g\ cm\ s^{-1}}$]')
    ax10.yaxis.set_tick_params(which='both', labelbottom=True)
    ax12.legend(fontsize=6,loc='lower right')


    fig.tight_layout()
    fig.subplots_adjust(hspace=0.13, wspace=0.04)

    fig.savefig('energy_momentum.pdf',format='pdf')
    fig.savefig('energy_momentum.png')
    plt.show()



main()











































