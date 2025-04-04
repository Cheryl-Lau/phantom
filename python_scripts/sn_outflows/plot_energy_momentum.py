
'''
Plot the evolution of energy and momentum contained within spheres of different radii
'''

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import glob 
import os 


Myr = 1E6*365*24*60*60 
pc = 3.086E+18

# Sync to shock arrival time 
t0_shock_ff = 0.023
t0_shock_cf = 0.172
t0_shock_scf = 0.091

# Time range 
tmin = 0.0
tmax = 0.37



def plot_line(ax,x,y,colour,labelstr):

    x_sorted = sorted(x)
    y_sorted = [ydum for _, ydum in sorted(zip(x,y))]
    x_sorted.pop(0)
    y_sorted.pop(0)
    ax.plot(x_sorted,y_sorted,color=colour,linewidth=1,label=labelstr)    
    ax.set_yscale('log')
    ax.set_xlim([tmin,tmax])
#    ax.set_xticklabels([])
#    ax.set_yticklabels([])
    ax.grid('on', linestyle='--')
#    ax.tick_params(axis='both', which='both', bottom=True, top=True, left=True, right=True, labelbottom=True, labelleft=True)
    
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




def read_rawfiles(filepath_energ,filepath_momen,chnl_label,df_energ,df_momen):

    # Energy files 
    for filename in glob.glob(filepath_energ):  # over_time 
        time = np.loadtxt(filename,skiprows=1,max_rows=1,unpack=True)
        radius, ke, te, tote = np.loadtxt(filename,skiprows=3,unpack=True)
        df_add = pd.DataFrame({'chnl':chnl_label,'time':time/Myr,'radius':radius/pc,'ekin':ke,'etherm':te,'etot':tote})
        df_energ = pd.concat([df_energ,df_add],ignore_index=True)

    # Momentum files 
    for filename in glob.glob(filepath_momen):  # over_time 
        with open(filename, 'r') as file:
            lines = file.readlines()
            time = float(lines[1].strip())
            mom_cioffi = float(lines[3].strip())
            mom_kimmcen = float(lines[4].strip())
        radius, mom = np.loadtxt(filename,skiprows=6,unpack=True)
        df_add = pd.DataFrame({'chnl':chnl_label,'time':time/Myr,'radius':radius/pc,'momen':mom,'momen_cioffi':mom_cioffi,'momen_kimmcen':mom_kimmcen})
        df_momen = pd.concat([df_momen,df_add],ignore_index=True)

    # Store spheres radii info 
    global nrad, radii
    nrad = len(radius)
    radii = radius/pc

    return df_energ, df_momen



def combine_frames():

    df_energ = pd.DataFrame(columns=('chnl', 'time', 'radius', 'ekin', 'etherm', 'etot'))
    df_momen = pd.DataFrame(columns=('chnl', 'time', 'radius', 'momen', 'momen_cioffi', 'momen_kimmcen'))

    # Free-field 
    filepath_energ = 'free_field/energy_insphere_freefield_*.dat'
    filepath_momen = 'free_field/momentum_insphere_freefield_*.dat'
    df_energ, df_momen = read_rawfiles(filepath_energ, filepath_momen, -1, df_energ, df_momen)

    # Confined 
    filepath_energ = 'confined/energy_insphere_confined_*.dat'
    filepath_momen = 'confined/momentum_insphere_confined_*.dat'
    df_energ, df_momen = read_rawfiles(filepath_energ, filepath_momen, 0, df_energ, df_momen)

    # Semi-confined 1chnl
    filepath_energ = 'semi_confined/vary_chnlsize/chnl005/energy_insphere_semiconf_1chnl005_*.dat'
    filepath_momen = 'semi_confined/vary_chnlsize/chnl005/momentum_insphere_semiconf_1chnl005_*.dat'
    df_energ, df_momen = read_rawfiles(filepath_energ, filepath_momen, 1, df_energ, df_momen)

    # Semi-confined 4chnl 
    filepath_energ = 'semi_confined/vary_nchnl/4chnl/energy_insphere_semiconf_4chnl005_*.dat'
    filepath_momen = 'semi_confined/vary_nchnl/4chnl/momentum_insphere_semiconf_4chnl005_*.dat'
    df_energ, df_momen = read_rawfiles(filepath_energ, filepath_momen, 4, df_energ, df_momen)

    print(df_energ)
    print(df_momen)

    return df_energ, df_momen


def plot_energ(ax,rad_pc,energ_type,df_energ): 

    # Free-field 
    iselect = np.where(np.logical_and(df_energ['radius'] == rad_pc, df_energ['chnl']==-1))[0]
    time = df_energ['time'].iloc[iselect].to_numpy()
    energ = df_energ[energ_type].iloc[iselect].to_numpy()
    plot_line(ax,time-t0_shock_ff,energ,'blue','free-field')

    # Confined 
    iselect = np.where(np.logical_and(df_energ['radius'] == rad_pc, df_energ['chnl']==0))[0]
    time = df_energ['time'].iloc[iselect].to_numpy()
    energ = df_energ[energ_type].iloc[iselect].to_numpy()
    plot_line(ax,time-t0_shock_cf,energ,'green','confined')

    # Semi-confined 1chnl
    iselect = np.where(np.logical_and(df_energ['radius'] == rad_pc, df_energ['chnl']==1))[0]
    time = df_energ['time'].iloc[iselect].to_numpy()
    energ = df_energ[energ_type].iloc[iselect].to_numpy()
    plot_line(ax,time-t0_shock_scf,energ,'red','semi-confined w/ 1 channel')

    # Semi-confined 4chnl
    iselect = np.where(np.logical_and(df_energ['radius'] == rad_pc, df_energ['chnl']==4))[0]
    time = df_energ['time'].iloc[iselect].to_numpy()
    energ = df_energ[energ_type].iloc[iselect].to_numpy()
    plot_dashed_line(ax,time-t0_shock_scf,energ,'red','semi-confined w/ 4 channels')

    return 


def plot_momen(ax,rad_pc,df_momen): 

    # Free-field 
    iselect = np.where(np.logical_and(df_momen['radius'] == rad_pc, df_momen['chnl']==-1))[0]
    time = df_momen['time'].iloc[iselect].to_numpy()
    momen = df_momen['momen'].iloc[iselect].to_numpy()
    momen_cioffi = df_momen['momen_cioffi'].iloc[iselect].to_numpy()
    momen_kimmcen = df_momen['momen_kimmcen'].iloc[iselect].to_numpy()
    plot_line(ax,time-t0_shock_ff,momen,'blue','free-field')
    plot_dotted_line(ax,time-t0_shock_ff,momen_cioffi,'black','Cioffi+ 1988')
    plot_dotted_line(ax,time-t0_shock_ff,momen_kimmcen,'grey','Kimm & Cen 2014')

    # Confined 
    iselect = np.where(np.logical_and(df_momen['radius'] == rad_pc, df_momen['chnl']==0))[0]
    time = df_momen['time'].iloc[iselect].to_numpy()
    momen = df_momen['momen'].iloc[iselect].to_numpy()
    plot_line(ax,time-t0_shock_cf,momen,'green','confined')

    # Semi-confined 1chnl
    iselect = np.where(np.logical_and(df_momen['radius'] == rad_pc, df_momen['chnl']==1))[0]
    time = df_momen['time'].iloc[iselect].to_numpy()
    momen = df_momen['momen'].iloc[iselect].to_numpy()
    plot_line(ax,time-t0_shock_scf,momen,'red','semi-confined w/ 1 channel')

    # Semi-confined 4chnl
    iselect = np.where(np.logical_and(df_momen['radius'] == rad_pc, df_momen['chnl']==4))[0]
    time = df_momen['time'].iloc[iselect].to_numpy()
    momen = df_momen['momen'].iloc[iselect].to_numpy()
    plot_dashed_line(ax,time-t0_shock_scf,momen,'red','semi-confined w/ 4 channels')

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

    return 



def main():

    df_energ, df_momen = combine_frames()

    fig = plt.figure(figsize=[10,10],dpi=100)

    # Total energy
    ax1 = fig.add_subplot(451)
    ax2 = fig.add_subplot(452)
    ax3 = fig.add_subplot(453)
    ax4 = fig.add_subplot(454)
    ax5 = fig.add_subplot(455)
    axes = [ax1,ax2,ax3,ax4,ax5]
    plot_energ_for_all_rad(axes,'etot',df_energ)
    for ax in axes:
        ax.set_ylim([1.5E49,5E51])
    ax1.set_ylabel('total energy [erg]')
    ax1.yaxis.set_tick_params(which='both', labelbottom=True)
    ax1.set_title('r = 15 pc')
    ax2.set_title('r = 30 pc')
    ax3.set_title('r = 45 pc')
    ax4.set_title('r = 60 pc')
    ax5.set_title('r = 75 pc')
    ax5.legend(fontsize=6.5,loc='lower right')

    # Kinetic energy 
    ax6 = fig.add_subplot(456)
    ax7 = fig.add_subplot(457)
    ax8 = fig.add_subplot(458)
    ax9 = fig.add_subplot(459)
    ax10 = fig.add_subplot(4,5,10)
    axes = [ax6,ax7,ax8,ax9,ax10]
    plot_energ_for_all_rad(axes,'ekin',df_energ)
    for ax in axes:
        ax.set_ylim([8E45,5E52])
    ax6.set_ylabel('kinetic energy [erg]')
    ax6.yaxis.set_tick_params(which='both', labelbottom=True)
    ax10.legend(fontsize=6.5,loc='lower right')

    # Thermal energy 
    ax11 = fig.add_subplot(4,5,11)
    ax12 = fig.add_subplot(4,5,12)
    ax13 = fig.add_subplot(4,5,13)
    ax14 = fig.add_subplot(4,5,14)
    ax15 = fig.add_subplot(4,5,15)
    axes = [ax11,ax12,ax13,ax14,ax15]
    plot_energ_for_all_rad(axes,'etherm',df_energ)
    for ax in axes:
        ax.set_ylim([8E48,5E51])
    ax11.set_ylabel('thermal energy [erg]')
    ax11.yaxis.set_tick_params(which='both', labelbottom=True)
    ax15.legend(fontsize=6.5,loc='lower right')

    # Momentum 
    ax16 = fig.add_subplot(4,5,16)
    ax17 = fig.add_subplot(4,5,17)
    ax18 = fig.add_subplot(4,5,18)
    ax19 = fig.add_subplot(4,5,19)
    ax20 = fig.add_subplot(4,5,20)
    axes = [ax16,ax17,ax18,ax19,ax20]
    plot_momen_for_all_rad(axes,df_momen)
    for ax in axes:
        ax.set_ylim([1E41,4E45])
    ax16.set_xlabel('time [Myr]')
    ax17.set_xlabel('time [Myr]')
    ax18.set_xlabel('time [Myr]')
    ax19.set_xlabel('time [Myr]')
    ax20.set_xlabel('time [Myr]')
    ax16.set_ylabel('radial momentum [$\mathrm{g\ cm\ s^{-1}}$]')
    ax16.yaxis.set_tick_params(which='both', labelbottom=True)
    ax20.legend(fontsize=6,loc='lower right')


    fig.tight_layout()
    fig.subplots_adjust(hspace=0.13, wspace=0.04)

    fig.savefig('energy_momentum_adiabatic.pdf',format='pdf')
    fig.savefig('energy_momentum_adiabatic.png')
    plt.show()



main()











































