#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import glob
import h5py
import pandas as pd 

udist_si = 3.085E16

for fname in sorted(glob.glob("snapshot010.hdf5")):
    file = h5py.File(fname,"r")

    box = np.array(file["/Header"].attrs["BoxSize"])
    box_center = 0.5 * box
    
    coords = np.array(file["/PartType0/Coordinates"])
    nfracH = np.array(file["/PartType0/NeutralFractionH"])

    x = coords[:,0]
    y = coords[:,1]
    z = coords[:,2]
    
    x_corr = coords[:,0]/udist_si - box_center[0]/udist_si
    y_corr = coords[:,1]/udist_si - box_center[1]/udist_si
    z_corr = coords[:,2]/udist_si - box_center[2]/udist_si
    dataset = pd.DataFrame({'x':x_corr, 'y':y_corr, 'z':z_corr, 'nH':nfracH}, columns=['x','y','z','nH'])
    dataset.to_csv('xyzf_cmi.txt',float_format='%0.4e', sep='\t')
    
    x_slice = []
    y_slice = []
    nH_slice = []
    for i in range(len(x)):
        
        if (z[i] > 0*udist_si+box_center[2] and z[i] < 0.2*udist_si+box_center[2]):
            x_slice.append(x[i]/udist_si)
            y_slice.append(y[i]/udist_si)
            nH_slice.append(nfracH[i])
            

    print('nH_slice',nH_slice)
    nH_scaled = nH_slice/max(nH_slice)
    
    plt.scatter(x_slice,y_slice,s=1,c=nH_scaled)
    

    for frac in nH_slice:
        if (frac < 0.):
            print('negative',frac)
    
    file.close()
            
    











