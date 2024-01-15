#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Reads nixyzhmf_* files to plot ionization radius over time
Compares results to Spitzer solution
'''

import numpy as np
import matplotlib.pyplot as plt 
import math
import glob


# Params in cgs units 
ncell = 5000
alpha = 2.7E-13
lumin = 1E49
rhozero = 5.21E-21
mass_p = 1.67262158E-24
cs_ion = 12.85E5
udist = 3.086E+18


def main():
    '''
    Loop over all nixyzhmf_* output files to plot evolution of 
    ionization front 
    '''
    time_ev = []
    ionrad_ev = []
    ionrad_spitz_ev = []
    
    for filename in glob.glob("nixyzhmf_0*.txt"):  # nixyzhmf_0*.txt
        print(filename)
        file = open(filename)
        for iline,line in enumerate(file):
            if (iline == 1):
                time_cgs = line 
                break
        time_cgs = float(time_cgs)

        nixyzhmf = np.loadtxt(filename,skiprows=3)
        ionrad = get_ionrad(nixyzhmf)
        
        # Calculate spitzer solution at time_cgs
        rad_strm = ((3.*lumin*mass_p**2)/(4.*np.pi*alpha*rhozero**2))**(1./3.)
        rad_spitz = rad_strm * (1. + 7./4.*(cs_ion*time_cgs/rad_strm))**(4./7.)
        
        print(ionrad/udist,rad_spitz/udist)
        
        # Convert to Myr,pc
        time_myr = time_cgs/(1E6*365*24*60*60) 
        ionrad_pc = ionrad/3.08567758128E+18
        rad_spitz_pc = rad_spitz/3.08567758128E+18
        
        if time_myr < 0.05:
            time_ev.append(time_myr)
            ionrad_ev.append(ionrad_pc)
            ionrad_spitz_ev.append(rad_spitz_pc)
    
    plt.plot(time_ev,ionrad_ev,'o',label='Phantom+CMI')
    plt.plot(time_ev,ionrad_spitz_ev,'o',label='Spitzer solution')
    plt.legend()
    plt.show()



def get_ionrad(nixyzhmf):
    '''
    Find ionization front from one output file 
    Input: 
        - nixyzhmf, node properties 
    Output:
        - ionrad, radius of ionization front 
    '''
    x_points = nixyzhmf[:,2]
    y_points = nixyzhmf[:,3]
    z_points = nixyzhmf[:,4]
    nH_points = nixyzhmf[:,7]
    
    r_points = np.sqrt(x_points**2+y_points**2+z_points**2)
    rmax = np.max(r_points) 
    
    # Sort onto 1-D grid of nH vs r 
    nHgrid = np.zeros(ncell)
    npgrid = np.zeros(ncell)

    rgrid = np.linspace(0,rmax,num=ncell)
    for ipoint in range(len(x_points)):
        r = r_points[ipoint]
        nH = nH_points[ipoint]
        icell = math.floor(r/rmax*ncell) - 1 
        nHgrid[icell] = nHgrid[icell] + nH
        npgrid[icell] = npgrid[icell] + 1
        
    nHgrid = np.divide(nHgrid,npgrid)

    # Interpolate to give full curve
    nans,nonzerofunc = nan_helper(nHgrid)
    nHgrid[nans]= np.interp(nonzerofunc(nans), nonzerofunc(~nans), nHgrid[~nans])
    
    # Find r between nH of 0.2 to 0.8 
    iclose02 = (np.abs(nHgrid - 0.2)).argmin()
    iclose08 = (np.abs(nHgrid - 0.8)).argmin()
    r02 = rgrid[iclose02]
    r08 = rgrid[iclose08]
    
    ionrad = (r02+r08)/2.
    ionrad = ionrad * udist  # convert to cgs units 
    
    return ionrad
    

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.
    https://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array
    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]
    
    
    
    




if __name__ == "__main__":
    main()
        
        
        
        
        
    
    












