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


for fname in sorted(glob.glob("snapshot010.hdf5")):
    file = h5py.File(fname,"r")

    file.visititems(printall)

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

    dataset = pd.DataFrame({'x':x_corr, 'y':y_corr, 'z':z_corr, 'rho':rho_cgs}, columns=['x','y','z','rho'])
    print(dataset)
    dataset.to_csv('xyzp_cmi.txt',float_format='%0.4e', sep='\t')

    file.close()
            
    











