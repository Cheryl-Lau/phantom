# -*- coding: utf-8 -*-

import numpy as np
from numpy.fft import fftn
from numpy import sqrt, zeros, conj, pi, arange, ones, convolve
import matplotlib.pyplot as plt 
import sys

Myr = 1e6*365*24*60*60
utime = 4.706E+14
t0_scf = 5.8810E-03*utime/Myr  # starting time of SN injection (exclude HII evolve time)
t0_ff = 1.4120E-03*utime/Myr


def read_dump_veldata(filename_vx,filename_vy,filename_vz):

    time_cgs = np.loadtxt(filename_vx,skiprows=1,max_rows=1)
    ix,iy,iz,x,y,z,vx = np.loadtxt(filename_vx,skiprows=2,unpack=True)
    ix,iy,iz,x,y,z,vy = np.loadtxt(filename_vy,skiprows=2,unpack=True)
    ix,iy,iz,x,y,z,vz = np.loadtxt(filename_vz,skiprows=2,unpack=True)

    nx = int(np.max(ix))
    ny = int(np.max(iy))
    nz = int(np.max(iz))

    if((nx != ny) or (ny != nz)):
        print(filename_vx,filename_vy,filename_vz)
        print(nx,ny,nz)
        print('Cube is not of the right size!!!!')
        exit()

    sizex = np.max(x)-np.min(x)
    sizey = np.max(y)-np.min(y)
    sizez = np.max(z)-np.min(z)

    vx_cube = zeros((nx,ny,nz))
    vy_cube = zeros((nx,ny,nz))
    vz_cube = zeros((nx,ny,nz))
    for iref in range(len(ix)):
        iix = int(ix[iref])-1
        iiy = int(iy[iref])-1
        iiz = int(iz[iref])-1
        vx_cube[iix,iiy,iiz] = vx[iref] 
        vy_cube[iix,iiy,iiz] = vy[iref]
        vz_cube[iix,iiy,iiz] = vz[iref]

    return time_cgs,vx_cube,vy_cube,vz_cube,sizex,sizey,sizez,nx,ny,nz
    

def Helmholtz_Hodge_decomposition(u,v,w,nx,ny,nz):
    '''
    Given three 3D velocity field, decompose into solenoidal and compressive parts
    source: https://github.com/shixun22/helmholtz/blob/master/helmholtz.py
    '''
    print('nx ny nz',nx,ny,nz)

    vx_f = np.fft.fftn(u)
    vy_f = np.fft.fftn(v)
    vz_f = np.fft.fftn(w)

    kx = np.fft.fftfreq(nx).reshape(nx,1,1)
    ky = np.fft.fftfreq(ny).reshape(ny,1)
    kz = np.fft.fftfreq(nz)

    k2 = kx**2 + ky**2 + kz**2
    k2[0,0,0] = 1. # to avoid inf. we do not care about the k=0 component

    div_Vf_f = (vx_f * kx +  vy_f * ky + vz_f * kz) #* 1j
    V_compressive_overk = div_Vf_f / k2
    V_compressive_x = np.fft.ifftn(V_compressive_overk * kx) #[:,np.newaxis,np.newaxis])
    V_compressive_y = np.fft.ifftn(V_compressive_overk * ky)
    V_compressive_z = np.fft.ifftn(V_compressive_overk * kz)

    V_solenoidal_x = u - V_compressive_x
    V_solenoidal_y = v - V_compressive_y
    V_solenoidal_z = w - V_compressive_z

    # check if the solenoidal part really divergence-free
    divVs = np.fft.ifftn((np.fft.fftn(V_solenoidal_x) * kx + np.fft.fftn(V_solenoidal_y) * ky + np.fft.fftn(V_solenoidal_z) * kz) * 1j * 2. * np.pi)

    print('div_solenoidal max:', abs(divVs).max())

    return V_compressive_x,V_compressive_y,V_compressive_z,V_solenoidal_x,V_solenoidal_y,V_solenoidal_z 



def fftfreq_doublenyquist(N):
    """
    fftfreq_ has only one negative Nyquist frequency -0.5
    fftfreq_doublenyquist has both -0.5 and 0.5
    """
    return np.arange(-(N/2), (N/2)+1) * 1. / N



def Helmholtz_Hodge_decomposition2(u,v,w,NX,NY,NZ):
    '''
    Given three 3D velocity field, decompose into solenoidal and compressive parts
    Slower but more accurate method 
    source: https://github.com/shixun22/helmholtz/blob/master/helmholtz_real.py
    '''

    V = np.vstack([u,v,w])

    vx_f = np.fft.fftn(u)
    vy_f = np.fft.fftn(v)
    vz_f = np.fft.fftn(w)
    x = np.arange(NX)
    y = np.arange(NY)
    z = np.arange(NZ)
    kx = fftfreq_doublenyquist(NX)
    ky = fftfreq_doublenyquist(NY)
    kz = fftfreq_doublenyquist(NZ)
    k2 = kx[:,np.newaxis,np.newaxis]**2 + ky[np.newaxis,:,np.newaxis]**2 + kz[np.newaxis,np.newaxis,:]**2
    k2[np.where(abs(k2)==0.)] = 1.

    V_k = np.array([[[[(V[xyz] * np.exp(-2.* np.pi * 1j * (kx[ikx] * x[:,np.newaxis,np.newaxis] + ky[iky] * y[np.newaxis,:,np.newaxis] + kz[ikz] * z[np.newaxis,np.newaxis,:]))).sum() for ikz in range(kz.size)] for iky in range(ky.size)] for ikx in range(kx.size)] for xyz in range(3)])

    if (NX%2)==0:
        V_k[:,0,:,:] = V_k[:,0,:,:]/2.
        V_k[:,-1,:,:] = V_k[:,-1,:,:]/2.

    if (NY%2)==0:
        V_k[:,:,0,:] = V_k[:,:,0,:]/2.
        V_k[:,:,-1,:] = V_k[:,:,-1,:]/2.

    if (NZ%2)==0:
        V_k[:,:,:,0] = V_k[:,:,:,0]/2.
        V_k[:,:,:,-1] = V_k[:,:,:,-1]/2.

    ## recontruct V from inverse Fourier transformation
    #V_reconstruct = np.array([[[[(V_k[xyz] * np.exp(2.* np.pi * 1j * (kx[:,np.newaxis,np.newaxis] * x[ix] + ky[np.newaxis,:,np.newaxis] * y[iy] + kz[np.newaxis,np.newaxis,:] * z[iz]))).sum() / NX / NY / NZ for iz in range(z.size)] for iy in range(y.size)] for ix in range(x.size)] for xyz in range(3)])

    div_V_k = (V_k[0] * kx[:,np.newaxis,np.newaxis] + V_k[1] * ky[np.newaxis,:,np.newaxis] + V_k[2] * kz[np.newaxis,np.newaxis,:]) / k2

    V_solenoidal_k = np.array([V_k[0] - div_V_k * kx[:,np.newaxis,np.newaxis], V_k[1] - div_V_k * ky[np.newaxis,:,np.newaxis], V_k[2] - div_V_k * kz[np.newaxis,np.newaxis,:]])

    V_solenoidal = np.array([[[[(V_solenoidal_k[xyz] * np.exp(2.* np.pi * 1j * (kx[:,np.newaxis,np.newaxis] * x[ix] + ky[np.newaxis,:,np.newaxis] * y[iy] + kz[np.newaxis,np.newaxis,:] * z[iz]))).sum() / NX / NY / NZ for iz in range(z.size)] for iy in range(y.size)] for ix in range(x.size)] for xyz in range(3)])


    # check if the solenoidal part is real 
    print('max of imaginary part of V_solenoidal:', abs(V_solenoidal.imag).max())

    # check if the solenoidal part really divergence-free
    divVs_r = np.array([[[((V_solenoidal_k[0] * kx[:,np.newaxis,np.newaxis] + V_solenoidal_k[1] * ky[np.newaxis,:,np.newaxis] + V_solenoidal_k[2] * kz[np.newaxis,np.newaxis,:]) / k2 * 2. * np.pi * 1j * np.exp(2.* np.pi * 1j * (kx[:,np.newaxis,np.newaxis] * x[ix] + ky[np.newaxis,:,np.newaxis] * y[iy] + kz[np.newaxis,np.newaxis,:] * z[iz]))).sum() / NX / NY / NZ for iz in range(z.size)] for iy in range(y.size)] for ix in range(x.size)])

    print('max of divergence(V_solenoidal):', abs(divVs_r).max())


    V_solenoidal_x = V_solenoidal[0]
    V_solenoidal_y = V_solenoidal[1]
    V_solenoidal_z = V_solenoidal[2]
    V_compressive_x = u - V_solenoidal_x
    V_compressive_y = v - V_solenoidal_y
    V_compressive_z = w - V_solenoidal_z

    return V_compressive_x,V_compressive_y,V_compressive_z,V_solenoidal_x,V_solenoidal_y,V_solenoidal_z


def compute_tke_spectrum(u,v,w,lx,ly,lz,smooth):
    """
    Given a velocity field u, v, w, this function computes the kinetic energy
    spectrum of that velocity field in spectral space. This procedure consists of the 
    following steps:
    1. Compute the spectral representation of u, v, and w using a fast Fourier transform.
    This returns uf, vf, and wf (the f stands for Fourier)
    2. Compute the point-wise kinetic energy Ef (kx, ky, kz) = 1/2 * (uf, vf, wf)* conjugate(uf, vf, wf)
    3. For every wave number triplet (kx, ky, kz) we have a corresponding spectral kinetic energy 
    Ef(kx, ky, kz). To extract a one dimensional spectrum, E(k), we integrate Ef(kx,ky,kz) over
    the surface of a sphere of radius k = sqrt(kx^2 + ky^2 + kz^2). In other words
    E(k) = sum( E(kx,ky,kz), for all (kx,ky,kz) such that k = sqrt(kx^2 + ky^2 + kz^2) ).

    Parameters:
    -----------  
    u: 3D array
    The x-velocity component.
    v: 3D array
    The y-velocity component.
    w: 3D array
    The z-velocity component.    
    lx: float
    The domain size in the x-direction.
    ly: float
    The domain size in the y-direction.
    lz: float
    The domain size in the z-direction.
    smooth: boolean
    A boolean to smooth the computed spectrum for nice visualization.

    source: https://turbulence.utah.edu/code/
    """
    nx = len(u[:,0,0])
    ny = len(v[0,:,0])
    nz = len(w[0,0,:])

    nt= nx*ny*nz
    n = nx #int(np.round(np.power(nt,1.0/3.0)))

    uh = fftn(u)/nt
    vh = fftn(v)/nt
    wh = fftn(w)/nt

    tkeh = zeros((nx,ny,nz))
    tkeh = 0.5*(uh*conj(uh) + vh*conj(vh) + wh*conj(wh)).real

    k0x = 2.0*pi/lx
    k0y = 2.0*pi/ly
    k0z = 2.0*pi/lz

    knorm = (k0x + k0y + k0z)/3.0

    kxmax = nx/2
    kymax = ny/2
    kzmax = nz/2

#    wave_numbers = knorm*arange(0,n)
    wave_numbers = knorm*np.logspace(0,np.log10(n),n)

    tke_spectrum = zeros(len(wave_numbers))

    for kx in range(nx):
        rkx = kx
        if (kx > kxmax):
            rkx = rkx - (nx)
        for ky in range(ny):
            rky = ky
            if (ky > kymax):
                rky = rky - (ny)
            for kz in range(nz):        
                rkz = kz
                if (kz > kzmax):
                    rkz = rkz - (nz)
                rk = sqrt(rkx*rkx + rky*rky + rkz*rkz)
                k = int(np.round(rk))
                tke_spectrum[k] = tke_spectrum[k] + tkeh[kx,ky,kz]

    tke_spectrum = tke_spectrum/knorm
    #  tke_spectrum = tke_spectrum[1:]
    #  wave_numbers = wave_numbers[1:]
    if (smooth == True):
        tkespecsmooth = movingaverage(tke_spectrum, 5) #smooth the spectrum
        tkespecsmooth[0:4] = tke_spectrum[0:4] # get the first 4 values from the original data
        tke_spectrum = tkespecsmooth

    knyquist = knorm*min(nx,ny,nz)/2 

    return knyquist,wave_numbers,tke_spectrum


def movingaverage(interval, window_size):
    window = ones(int(window_size))/float(window_size)

    return convolve(interval, window, 'same')


def plot_powerspec(ax,time,wavenumber,spectrum,linestylestr,colourstr,labelstr):

    ax.plot(wavenumber,spectrum,linestylestr.strip(),color=colourstr.strip(),label=labelstr.strip())
    ax.axvspan(5.6e0,2e1,alpha=0.3,color='lightgrey')
    ax.axvspan(2e-1,4e-1,alpha=0.3,color='lightgrey')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('wave number k [$\mathrm{pc^{-1}}$]')
    ax.set_ylabel('kinetic energy E(k) [erg]')
    ax.set_xlim([2.7e-1,9e0])
    ax.set_ylim([5e5,5e16])
    ax.legend(loc='upper right',fontsize=7)

    return 


def create_spec(ax,filename_vx,filename_vy,filename_vz,labelstr,colourstr,textloc_x,textloc_y):

    time_cgs,vx_cube,vy_cube,vz_cube,sizex,sizey,sizez,nx,ny,nz = read_dump_veldata(filename_vx,filename_vy,filename_vz)

    vx_compr,vy_compr,vz_compr,vx_solen,vy_solen,vz_solen = Helmholtz_Hodge_decomposition(vx_cube,vy_cube,vz_cube,nx,ny,nz)

    knyquist, wavenumber,tke_spectrum = compute_tke_spectrum(vx_cube,vy_cube,vz_cube,sizex,sizey,sizez,True)
    plot_powerspec(ax,time_cgs,wavenumber,tke_spectrum,':',colourstr,labelstr+' total $E_\mathrm{kin}$')

    knyquist, wavenumber,tke_spectrum = compute_tke_spectrum(vx_compr,vy_compr,vz_compr,sizex,sizey,sizez,True)
    plot_powerspec(ax,time_cgs,wavenumber[1:],tke_spectrum[1:],'-',colourstr,labelstr+' compressive mode')

    knyquist, wavenumber,tke_spectrum = compute_tke_spectrum(vx_solen,vy_solen,vz_solen,sizex,sizey,sizez,True)
    plot_powerspec(ax,time_cgs,wavenumber,tke_spectrum,'--',colourstr,labelstr+' solenoidal mode')


    if (labelstr == 'semi-confined'):
        time_Myr = time_cgs/Myr - t0_scf
    else:
        time_Myr = time_cgs/Myr - t0_ff
    ax.text(textloc_x,textloc_y,labelstr+': '+str(round(time_Myr,3))+' Myr',fontsize=7)

    return 


def main():

    fig1,(ax1,ax2,ax3,ax4) = plt.subplots(nrows=4,sharex=False,subplot_kw=dict(frameon=True),figsize=(5,12),dpi=200)


    # Shock front 15 pc: 
    ax1.text(3.4e-1,3e15,'Shock front at 15 pc')

    # Semi-confined case 
    filename_vx = 'semi_confined/vary_chnlsize/chnl005/velfield_x_semiconf_1chnl005_00250.dat'
    filename_vy = 'semi_confined/vary_chnlsize/chnl005/velfield_y_semiconf_1chnl005_00250.dat'
    filename_vz = 'semi_confined/vary_chnlsize/chnl005/velfield_z_semiconf_1chnl005_00250.dat'
    create_spec(ax1,filename_vx,filename_vy,filename_vz,'semi-confined','red',3.4e-1,7e14)

    # Free-field case 
    filename_vx = 'free_field/velfield_x_freefield_00003.dat'
    filename_vy = 'free_field/velfield_y_freefield_00003.dat'
    filename_vz = 'free_field/velfield_z_freefield_00003.dat'
    create_spec(ax1,filename_vx,filename_vy,filename_vz,'free-field','darkblue',3.4e-1,2e14)


    # Shock front 30 pc: 
    ax2.text(3.4e-1,3e15,'Shock front at 30 pc')

    # Semi-confined case 
    filename_vx = 'semi_confined/vary_chnlsize/chnl005/velfield_x_semiconf_1chnl005_00650.dat'
    filename_vy = 'semi_confined/vary_chnlsize/chnl005/velfield_y_semiconf_1chnl005_00650.dat'
    filename_vz = 'semi_confined/vary_chnlsize/chnl005/velfield_z_semiconf_1chnl005_00650.dat'
    create_spec(ax2,filename_vx,filename_vy,filename_vz,'semi-confined','red',3.4e-1,7e14)

    # Free-field case 
    filename_vx = 'free_field/velfield_x_freefield_00014.dat'
    filename_vy = 'free_field/velfield_y_freefield_00014.dat'
    filename_vz = 'free_field/velfield_z_freefield_00014.dat'
    create_spec(ax2,filename_vx,filename_vy,filename_vz,'free-field','darkblue',3.4e-1,2e14)


    # Shock front 45 pc: 
    ax3.text(3.4e-1,3e15,'Shock front at 45 pc')

    # Semi-confined case 
    filename_vx = 'semi_confined/vary_chnlsize/chnl005/velfield_x_semiconf_1chnl005_01500.dat'
    filename_vy = 'semi_confined/vary_chnlsize/chnl005/velfield_y_semiconf_1chnl005_01500.dat'
    filename_vz = 'semi_confined/vary_chnlsize/chnl005/velfield_z_semiconf_1chnl005_01500.dat'
    create_spec(ax3,filename_vx,filename_vy,filename_vz,'semi-confined','red',3.4e-1,7e14)

    # Free-field case 
    filename_vx = 'free_field/velfield_x_freefield_00038.dat'
    filename_vy = 'free_field/velfield_y_freefield_00038.dat'
    filename_vz = 'free_field/velfield_z_freefield_00038.dat'
    create_spec(ax3,filename_vx,filename_vy,filename_vz,'free-field','darkblue',3.4e-1,2e14)


    # Shock front at 60 pc: 
    ax4.text(3.4e-1,3e15,'Shock front at 60 pc')

    # Semi-confined case 
    filename_vx = 'semi_confined/vary_chnlsize/chnl005/velfield_x_semiconf_1chnl005_02300.dat' 
    filename_vy = 'semi_confined/vary_chnlsize/chnl005/velfield_y_semiconf_1chnl005_02300.dat'
    filename_vz = 'semi_confined/vary_chnlsize/chnl005/velfield_z_semiconf_1chnl005_02300.dat'
    create_spec(ax4,filename_vx,filename_vy,filename_vz,'semi-confined','red',3.4e-1,7e14)

    # Free-field case 
    filename_vx = 'free_field/velfield_x_freefield_00080.dat'
    filename_vy = 'free_field/velfield_y_freefield_00080.dat'
    filename_vz = 'free_field/velfield_z_freefield_00080.dat'
    create_spec(ax4,filename_vx,filename_vy,filename_vz,'free-field','darkblue',3.4e-1,2e14)


    fig1.tight_layout(pad=0.5)
    fig1.savefig('turb_spec_adiabatic.pdf',format='pdf')
    fig1.savefig('turb_spec_adiabatic.png',dpi=200)

    plt.show()



if __name__ == "__main__":
    main()










