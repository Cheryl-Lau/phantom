
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import pandas as pd 
import glob 
from pathlib import Path


utime = 4.706e14
udist = 3.086e+19
umass = 1.989e+33
unit_density = 6.768e-23
unit_velocity = 6.558e+3
unit_ergg = 4.301e+7
Myr = 1e6*365*24*60*60 
solarm = 1.989e+33



def plot_settings(ax,ylabel,ymin,ymax):

    ax.set_yscale('log')

    ax.set_xlim([-0.001,0.061])
    ax.set_xlabel(r'$t_\mathrm{SN}$ [Myr]')

    ax.set_ylim([ymin,ymax])
    ax.set_ylabel(ylabel)

    ax.grid(linestyle='--',color='darkgrey',linewidth=0.5)

    return 


def bin_array(A_arr,T_arr,Amin,Amax,N):

    logAmin = np.log10(Amin)
    logAmax = np.log10(Amax)
    A_binned = np.logspace(logAmin,logAmax,N)
    dA = (logAmax-logAmin)/N 

    T_binned = [0.0]*N
    nT_binned = [0]*N

    for A,T in zip(A_arr,T_arr):
        iAbin = int((np.log10(A) - logAmin + dA) /dA)-1
        if (iAbin >= N):
            iAbin = N-1
        elif (iAbin < 0):
            iAbin = 0
        T_binned[iAbin] += T
        nT_binned[iAbin] += 1 

    for i in range(len(T_binned)):
        if (nT_binned[i] > 0):
            T_binned[i] = T_binned[i]/nT_binned[i]
    
    return np.array(A_binned), np.array(T_binned)



def get_mesh(sdir,A,Amin,Amax):

    path = sdir+"/outflow_6pc/gasflow_shell_*"

    meshxfile_path = Path('outflow_meshes_6pc/'+sdir+'_'+A+'_mesh_x.dat')
    meshyfile_path = Path('outflow_meshes_6pc/'+sdir+'_'+A+'_mesh_y.dat')
    meshzfile_path = Path('outflow_meshes_6pc/'+sdir+'_'+A+'_mesh_z.dat')

    if meshxfile_path.exists() and meshyfile_path.exists() and meshyfile_path.exists() :

        x_mesh = np.loadtxt(meshxfile_path)
        y_mesh = np.loadtxt(meshyfile_path)
        z_mesh = np.loadtxt(meshzfile_path)

    else:

        nbin = int(1e4)
        ntimefile = len(glob.glob(path))
        x_mesh = np.zeros((ntimefile))
        y_mesh = np.zeros((nbin))
        z_mesh = np.zeros((nbin,ntimefile))

        for it,filename in enumerate(sorted(glob.glob(path))): 
            print(filename)
        
            time = np.loadtxt(filename,skiprows=1,max_rows=1)
            if (it==0):
                t_SN = time/Myr
            time_Myr = time/Myr - t_SN

            data = np.loadtxt(filename,skiprows=3)

            df = pd.DataFrame(data, columns=['ip','T','rho','radvel','radmomen','pressure','rampress'])


            # Take the *change* in properties 
            if (it==0):
                df0 = df 
            else:
                df['radvel'] = df['radvel'] - df0['radvel']
                df['radmomen'] = df['radmomen'] - df0['radmomen']
                df['pressure'] = df['pressure'] - df0['pressure']
                df['rampress'] = df['rampress'] - df0['rampress']


            # Only those which are outflowing and is increased
            iout = np.where((df['radvel'] > 0) & (df['radmomen'] > 0) & (df['pressure'] > 0) & (df['rampress'] > 0))[0]
            df_out = df.loc[iout]

            A_binned,T_binned = bin_array(df_out[A],df_out['T'],Amin,Amax,nbin)

            x_mesh[it] = time_Myr
            y_mesh[0:nbin] = A_binned[0:nbin]
            z_mesh[0:nbin,it] = T_binned[0:nbin]

        np.savetxt(meshxfile_path,x_mesh)
        np.savetxt(meshyfile_path,y_mesh)
        np.savetxt(meshzfile_path,z_mesh)
        print('Mesh written to file for ',A)

    return x_mesh, y_mesh, z_mesh



def plot_analyt_model(ax,strA,filepath,t_SN):

    time_cgs,vx_cgs,thermpr_cgs,rampr_cgs = np.loadtxt(filepath,unpack=True)
    time_Myr = time_cgs/Myr + t_SN

    if (strA == 'radvel'):
        ax.plot(time_Myr,vx_cgs,'--',color='black',label='analytical model')
    elif (strA == 'pressure'):
        ax.plot(time_Myr,thermpr_cgs,'--',color='black',label='analytical model')
    elif (strA == 'rampress'):
        ax.plot(time_Myr,rampr_cgs,'--',color='black',label='analytical model')

    ax.legend(loc='upper right')

    return 



'''
simdir = ['turbcloud_semiconf','turbff_smear','turbff_env','ff_smear','ff_env']
simlabel = ['Semi-confined SN',r'Turb FF - $\rho_\mathrm{smear}$',r'Turb FF - $\rho_\mathrm{env}$',r'FF - $\rho_\mathrm{smear}$',r'FF - $\rho_\mathrm{env}$']
imfilename = ['outflow_cloud.png','outflow_tffsmear.png','outflow_tffenv.png','outflow_ffsmear.png','outflow_ffenv.png']
analytmodelfile = ['semi_confined_model.txt','free_field_smear_model.txt','free_field_env_model.txt','free_field_smear_model.txt','free_field_env_model.txt']
t_SN = [0.003,0.003,0.002,0.009,0]  # Myr 
'''
simdir = ['turbcloud_semiconf','turbff_smear']
simlabel = ['Semi-confined SN',r'Turb FF - $\rho_\mathrm{smear}$']
imfilename = ['outflow_cloud.png','outflow_tffsmear.png']
analytmodelfile = ['semi_confined_model.txt','free_field_smear_model.txt']
t_SN = [0.003,0.003]  # Myr


for sdir, modelfile, label, imfile, tSN in zip(simdir,analytmodelfile,simlabel,imfilename,t_SN):

    modelfile = 'analytical_model/at_6pc/'+modelfile

    fig, axes = plt.subplots(nrows=2, ncols=2, subplot_kw=dict(frameon=True), figsize=(10,8))
    ax_vr = axes[0,0]
    ax_mu = axes[1,0]
    ax_ptherm = axes[0,1]
    ax_pram = axes[1,1]

    vr_min = 5e4
    vr_max = 5e9
    x_mesh, y_mesh, z_mesh = get_mesh(sdir,'radvel',vr_min,vr_max)
    pcm = ax_vr.pcolormesh(x_mesh, y_mesh, np.log10(z_mesh), shading='nearest', cmap='viridis',vmin=1,vmax=5)
    plot_analyt_model(ax_vr,'radvel',modelfile,tSN)
    plot_settings(ax_vr,r'$v_r\ \mathrm{[cm\ s^{-1}]}$',vr_min,vr_max)

    mu_min = 5e30
    mu_max = 5e42
    x_mesh, y_mesh, z_mesh = get_mesh(sdir,'radmomen',mu_min,mu_max)
    pcm = ax_mu.pcolormesh(x_mesh, y_mesh, np.log10(z_mesh), shading='nearest', cmap='viridis',vmin=1,vmax=5)
    plot_settings(ax_mu,r'$\mu_r\ \mathrm{[g\ cm\ s^{-1}]}$',mu_min,mu_max)

    ptherm_min = 5e-18
    ptherm_max = 5e0
    x_mesh, y_mesh, z_mesh = get_mesh(sdir,'pressure',ptherm_min,ptherm_max)
    pcm = ax_ptherm.pcolormesh(x_mesh, y_mesh, np.log10(z_mesh), shading='nearest', cmap='viridis',vmin=1,vmax=5)
    plot_analyt_model(ax_ptherm,'pressure',modelfile,tSN)
    plot_settings(ax_ptherm,r'$p_\mathrm{therm}\ \mathrm{[g\ cm^{-1}\ s^{-2}]}$',ptherm_min,ptherm_max)

    pram_min = 5e-18
    pram_max = 5e0
    x_mesh, y_mesh, z_mesh = get_mesh(sdir,'rampress',pram_min,pram_max)
    pcm = ax_pram.pcolormesh(x_mesh, y_mesh, np.log10(z_mesh), shading='nearest', cmap='viridis',vmin=1,vmax=5)
    plot_analyt_model(ax_pram,'rampress',modelfile,tSN)
    plot_settings(ax_pram,r'$p_\mathrm{ram}\ \mathrm{[g\ cm^{-1}\ s^{-2}]}$',pram_min,pram_max)

    fig.suptitle(label,fontsize=15)
    fig.tight_layout()

    cbar = fig.colorbar(pcm, ax=(ax_vr,ax_mu,ax_ptherm,ax_pram), cmap='viridis', location='right', shrink=0.8,pad=0.03)
    cbar.set_label('log(T) [K]',rotation=270,labelpad=15)


    plt.savefig(imfile,dpi=200)

plt.show()














