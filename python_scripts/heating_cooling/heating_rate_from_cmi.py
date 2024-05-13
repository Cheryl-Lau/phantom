#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 18:42:28 2023

@author: okitacscl1
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
import glob


get_mean_gamma = False


times = []
gammas = []

for f in sorted(glob.glob("snapshot010.hdf5")):

    file = h5py.File(f, "r")
    
    datasetnames = file.keys()
    
    print(file['/Units'].attrs.keys())


    udist = file['/Units'].attrs['Unit length in cgs (U_L)']
    umass = file['/Units'].attrs['Unit mass in cgs (U_M)']
    utime = file['/Units'].attrs['Unit time in cgs (U_t)']

    unit_energ = umass*(udist/utime)**2
    
    print('unit_energ',unit_energ)
    
    unit_gamma = unit_energ/utime



    times.append(file["/Header"].attrs["Time"])

    coords = np.array(file["/PartType0/Coordinates"])
    neutfracH = np.array(file["/PartType0/NeutralFractionH"])

    gamma = np.array(file["/PartType0/HeatingRateH"])

    if get_mean_gamma == True:
        gamma_ion = []
        for i in range(len(neutfracH)):
            if neutfracH[i] > 0.1 and neutfracH[i] < 0.6:
                gamma_ion.append(gamma[i]*unit_gamma)  # to cgs
        
        gamma_mean = np.mean(gamma_ion)
        print(gamma_mean)
        
    else:
        ncell = 100
        gamma_cells = [0] * ncell
        num_cells = [0] * ncell
        nH_cells = np.linspace(0,1,ncell)
        
        for ip in range(len(neutfracH)):
            nH = neutfracH[ip]
            icell = np.abs(nH_cells - nH).argmin()
            gamma_cells[icell] += gamma[ip]*unit_gamma * neutfracH[ip]  # to cgs
            num_cells[icell] += 1
            
        for icell in range(ncell):
            if num_cells[icell] >= 1:
                gamma_cells[icell] = gamma_cells[icell] / num_cells[icell]
            else:
                gamma_cells[icell] = 0
        
        plt.plot(nH_cells,gamma_cells)
        plt.xlabel('nH')
        plt.ylabel('gamma [erg/s]')
        plt.yscale('log')
        plt.show()

    file.close()
    
    
