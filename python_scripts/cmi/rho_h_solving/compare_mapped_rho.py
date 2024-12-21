#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import glob
import h5py
import pandas as pd 

def printall(name, obj):
    print(name, dict(obj.attrs))


udist_si = 3.085E16
pmass = 1e-1*1.9891E+33
mass_proton_cgs = 1.67262158E-24
mass_proton_si = mass_proton_cgs*1E-3 


def read_snapshot(filename):

    file = h5py.File(filename,"r")

    box = np.array(file["/Header"].attrs["BoxSize"])
    box_center = 0.5 * box
    
    coords = np.array(file["/PartType0/Coordinates"])
    nfracH = np.array(file["/PartType0/NeutralFractionH"])
    numden = np.array(file["/PartType0/NumberDensity"])

    x = coords[:,0]
    y = coords[:,1]
    z = coords[:,2]
    
    x_corr = coords[:,0]/udist_si - box_center[0]/udist_si
    y_corr = coords[:,1]/udist_si - box_center[1]/udist_si
    z_corr = coords[:,2]/udist_si - box_center[2]/udist_si

    rho_si = numden*mass_proton_si
    rho_cgs = rho_si*1E-3

    df_snapshot = pd.DataFrame({'x':x_corr, 'y':y_corr, 'z':z_corr, 'rho':rho_cgs}, columns=['x','y','z','rho'])
    file.close()

    return df_snapshot  


def plot_rho_distri(ax,df_snapshot,minrho,maxrho,colourstr,labelstr): 

    rho = df_snapshot['rho'].to_numpy()
    ax.hist(rho,bins=np.logspace(np.log10(minrho),np.log10(maxrho),1000),histtype=u'step',color=colourstr,label=labelstr)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'density [$\mathrm{g\ cm^{-3}}$]')
    ax.set_ylabel('number of cells')
    ax.legend(loc='upper left')

    return 


def plot_h_distri(ax,filename,minh,maxh,colourstr,labelstr):

    h_vals = np.loadtxt(filename)
    ax.hist(h_vals,bins=np.logspace(np.log10(minh),np.log10(maxh),1000),histtype=u'step',color=colourstr,label=labelstr)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('smoothing length [pc]')
    ax.set_ylabel('number of pseudo-particles')
    ax.legend(loc='upper left')

    return



fig = plt.figure(figsize=(15,10))

ax1 = fig.add_subplot(231)
df_esth = read_snapshot('est_h/unifdis/snapshot010.hdf5')
df_solveh = read_snapshot('solve_h/unifdis/snapshot010.hdf5')
minrho = 3e-20
maxrho = 3e-19
plot_rho_distri(ax1,df_solveh,minrho,maxrho,'blue','solved iteratively')
plot_rho_distri(ax1,df_esth,minrho,maxrho,'cyan','estimated with node size')
ax1.set_title('uniform medium')
ax1.set_ylim([1E0,2E3])

ax2 = fig.add_subplot(232)
df_esth = read_snapshot('est_h/clumpy/snapshot010.hdf5')
df_solveh = read_snapshot('solve_h/clumpy/snapshot010.hdf5')
minrho = 2e-29
maxrho = 2e-17
plot_rho_distri(ax2,df_solveh,minrho,maxrho,'blue','solved iteratively')
plot_rho_distri(ax2,df_esth,minrho,maxrho,'cyan','estimated with node size')
ax2.set_title('clumpy medium')
ax2.set_ylim([1E0,2E3])

ax3 = fig.add_subplot(233)
df_esth = read_snapshot('est_h/very_clumpy/snapshot010.hdf5')
df_solveh = read_snapshot('solve_h/very_clumpy/snapshot010.hdf5')
minrho = 2e-30
maxrho = 2e-18
plot_rho_distri(ax3,df_solveh,minrho,maxrho,'blue','solved iteratively')
plot_rho_distri(ax3,df_esth,minrho,maxrho,'cyan','estimated with node size')
ax3.set_title('very clumpy medium')
ax3.set_ylim([1E0,6E2])

ax4 = fig.add_subplot(234)
minh = 1E-2
maxh = 0.15
plot_h_distri(ax4,'est_h/unifdis/smoothing_lengths.dat',minh,maxh,'cyan','estimated with node size')
plot_h_distri(ax4,'solve_h/unifdis/smoothing_lengths.dat',minh,maxh,'blue','solved iteratively')
ax4.set_ylim([1E0,2E3])

ax5 = fig.add_subplot(235)
minh = 4E-3
maxh = 0.25
plot_h_distri(ax5,'est_h/clumpy/smoothing_lengths.dat',minh,maxh,'cyan','estimated with node size')
plot_h_distri(ax5,'solve_h/clumpy/smoothing_lengths.dat',minh,maxh,'blue','solved iteratively')
ax5.set_ylim([1E0,2E3])

ax6 = fig.add_subplot(236)
minh = 0.01
maxh = 0.65
plot_h_distri(ax6,'est_h/very_clumpy/smoothing_lengths.dat',minh,maxh,'cyan','estimated with node size')
plot_h_distri(ax6,'solve_h/very_clumpy/smoothing_lengths.dat',minh,maxh,'blue','solved iteratively')
ax6.set_ylim([1E0,2E3])



plt.savefig('hmethod_compare_rhodistri.png',dpi=200)
plt.show()











