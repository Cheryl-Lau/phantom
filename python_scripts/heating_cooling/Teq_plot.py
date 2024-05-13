#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot Teqs wrt rho and gamma 
"""

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

plot_3D = True
get_slice = False 
plot_baro_eos = False 

data = np.loadtxt('rho_gamma_Teqs.txt',skiprows=1)

rho = data[:,0]
gamma = data[:,1]
Teq1 = data[:,3]
Teq2 = data[:,4]
Teq3 = data[:,5] 

# stores all ind points 
rho_all = []
gamma_all = []
Teq_all = []

for i in range(len(rho)):
    if (Teq1[i] != 0.):
        rho_all.append(rho[i])
        gamma_all.append(gamma[i])
        Teq_all.append(Teq1[i])
        
for i in range(len(rho)):
    if (Teq2[i] != 0.):
        rho_all.append(rho[i])
        gamma_all.append(gamma[i])
        Teq_all.append(Teq2[i])

for i in range(len(rho)):
    if (Teq3[i] != 0.):
        rho_all.append(rho[i])
        gamma_all.append(gamma[i])
        Teq_all.append(Teq3[i])
        

class MyAxes3D(axes3d.Axes3D):

    def __init__(self, baseObject, sides_to_draw):
        self.__class__ = type(baseObject.__class__.__name__,
                              (self.__class__, baseObject.__class__),
                              {})
        self.__dict__ = baseObject.__dict__
        self.sides_to_draw = list(sides_to_draw)
        self.mouse_init()

    def set_some_features_visibility(self, visible):
        for t in self.w_zaxis.get_ticklines() + self.w_zaxis.get_ticklabels():
            t.set_visible(visible)
        self.w_zaxis.line.set_visible(visible)
        self.w_zaxis.pane.set_visible(visible)
        self.w_zaxis.label.set_visible(visible)

    def draw(self, renderer):
        # set visibility of some features False 
        self.set_some_features_visibility(False)
        # draw the axes
        super(MyAxes3D, self).draw(renderer)
        # set visibility of some features True. 
        # This could be adapted to set your features to desired visibility, 
        # e.g. storing the previous values and restoring the values
        self.set_some_features_visibility(True)

        zaxis = self.zaxis
        draw_grid_old = zaxis.axes._draw_grid
        # disable draw grid
        zaxis.axes._draw_grid = False

        tmp_planes = zaxis._PLANES

        if 'l' in self.sides_to_draw :
            # draw zaxis on the left side
            zaxis._PLANES = (tmp_planes[2], tmp_planes[3],
                             tmp_planes[0], tmp_planes[1],
                             tmp_planes[4], tmp_planes[5])
            zaxis.draw(renderer)
        if 'r' in self.sides_to_draw :
            # draw zaxis on the right side
            zaxis._PLANES = (tmp_planes[3], tmp_planes[2], 
                             tmp_planes[1], tmp_planes[0], 
                             tmp_planes[4], tmp_planes[5])
            zaxis.draw(renderer)

        zaxis._PLANES = tmp_planes

        # disable draw grid
        zaxis.axes._draw_grid = draw_grid_old


if plot_3D == True:
    # 3D plot of all points 
    fig = plt.figure(figsize=[18,10],dpi=200)
    fig.subplots_adjust(wspace=0)
    
    ax1 = fig.add_subplot(122,projection='3d')
    
    ax1.set_ylabel(r'$ \mathrm{ log(\rho) \ [g \ cm^{-3}] }$')
    ax1.set_xlabel(r'$ \mathrm{ log(\Gamma) \ [erg \ s^{-1}] }$')
    ax1.invert_yaxis()
    ax1.zaxis.set_rotate_label(True)
    ax1.set_zlabel(r'$ \mathrm{ log(T_{eq}) \ [K] }$', rotation=90)
#    ax1.view_init(elev=10., azim=150.)
    ax1.view_init(elev=10., azim=330.)
    
#    ax1 = fig.add_axes(MyAxes3D(ax1, 'r'))
#    ax1.zaxis.set_ticks_position('top')
#    ax1.zaxis.set_label_position('top')
    
    pcm = ax1.scatter(np.log10(gamma_all), np.log10(rho_all), np.log10(Teq_all),c=np.log10(Teq_all), cmap=cm.bwr)
    
    # Get rid of the panes
    ax1.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax1.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax1.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        
    ax2 = fig.add_subplot(121,projection='3d')

    ax2.set_ylabel(r'$ \mathrm{ log(\rho) \ [g \ cm^{-3}] }$')
    ax2.set_xlabel(r'$ \mathrm{ log(\Gamma) \ [erg \ s^{-1}] }$')
    ax2.invert_yaxis()
    ax2.zaxis.set_rotate_label(False)
    ax2.set_zlabel(r'$ \mathrm{ log(T_{eq}) \ [K] }$', rotation=90)
#    ax2.view_init(elev=10., azim=60.)
    ax2.view_init(elev=10., azim=240.)
    
    pcm = ax2.scatter(np.log10(gamma_all), np.log10(rho_all), np.log10(Teq_all),c=np.log10(Teq_all), cmap=cm.bwr)
    
    ax2.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax2.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax2.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    
    cbar = fig.colorbar(pcm, ax=[ax1,ax2], location='right', shrink=0.5)
    cbar.set_label('log(T) [K]', rotation=270, labelpad=20)
    
    
    plt.savefig("Teq_solutions.png")

    


R_cgs = 8.31E7 

def barotrope(rho_cgs):
    
    if rho_cgs <= 5.5E-19:
        polyindex = 0.75
        pressure = 1.66E4 * rho_cgs**polyindex
    elif rho_cgs <= 5.5E-15:
        polyindex = 1.0
        pressure = 6.5E8 * rho_cgs**polyindex
    elif rho_cgs <= 2E-13:
        polyindex = 1.4
        pressure = 3E14 * rho_cgs**polyindex
    else:
        polyindex = 1.0
        pressure = 2.5E9 * rho_cgs**polyindex
    
    return pressure 


if get_slice == True:

    # Extract points for one particular gamma 
    gamma_extract = 1E-26
    tol_gamma = 1E-27
    rho_fixgamma = []
    Teq_fixgamma = []
    for i in range(len(rho_all)):
        if gamma_all[i] > (gamma_extract-tol_gamma) and gamma_all[i] < (gamma_extract+tol_gamma):
            rho_fixgamma.append(rho_all[i])
            Teq_fixgamma.append(Teq_all[i])
    fig2 = plt.figure(figsize=[7,5])
    ax2 = fig2.add_subplot(111)
    ax2.plot(rho_fixgamma, Teq_fixgamma,markersize=1,color='red',label='Heating & Cooling')
    ax2.set_xscale('log')
#    ax2.set_yscale('log')

    ax2.set_xlabel(r'$ \rho \ \mathrm{[g \ cm^{-3}]} $')
    ax2.set_ylabel(r'$T_{eq} \ \mathrm{[K]} $')
    
    # Barotropic EOS for comparison 
    if plot_baro_eos == True:
        rho_baro = np.logspace(-24,-11,100)
        pressure_all = []
        for rho in rho_baro:
            pressure = barotrope(rho)
            pressure_all.append(pressure)
    
        temperature_all = np.array(pressure_all)/(np.array(rho_baro)*R_cgs)
    
        ax2.plot(rho_baro,temperature_all,markersize=1,color='black',label='Barotropic EOS')
        
        
    ax2.set_xscale('log')
    
    if plot_baro_eos == False:
        ax2.set_yscale('log')
    else:
        ax2.set_ylim([0,200])
    
    ax2.legend()
    
    fig2.savefig('EOS_compare.png')


















