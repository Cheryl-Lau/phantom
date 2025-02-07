# -*- coding: utf-8 -*-
"""
Plots the evolution of ionization front radius 
compares on-tree to on-parts 
"""
import numpy as np
import matplotlib.pyplot as plt 
import glob
import os
import math

plot_tree = True
plot_parts = False
show_spitzer = True
show_hosoinu = True


# Params in cgs units 
ncell = 5000
alpha = 2.7E-13
lumin = 1E49
rhozero = 5.21E-21
mass_p = 1.67262158E-24
cs_ion = 12.85E5
udist = 3.086E+18


def main(): 
    
    fig, (ax1,ax2) = plt.subplots(nrows=2, sharex=True, subplot_kw=dict(frameon=True), gridspec_kw={'height_ratios': [4, 1]}, figsize=(7,6), dpi=300)
    plt.subplots_adjust(hspace=.0)
    
    if plot_tree == True:

        file_exist = os.path.isfile('./time_rad_spitz_tree.txt')
        if file_exist == True:
            # read existing file
            time_tree,ionrad_tree,spitzrad_tree,hosoinurad_tree = np.loadtxt('time_rad_spitz_tree.txt',unpack=True)
        else:
            # compute ion-rad
            time_tree,ionrad_tree,spitzrad_tree,hosoinurad_tree = extract_data(True)
            # save to file 
            data = np.column_stack([time_tree,ionrad_tree,spitzrad_tree,hosoinurad_tree])
            np.savetxt('time_rad_spitz_tree.txt',data)
            
        # Plot
        if show_spitzer == True and plot_parts == False:
            ax1.plot(time_tree,spitzrad_tree,color='black',label='Spitzer solution')
            errtree_spitzer = (ionrad_tree-spitzrad_tree)/spitzrad_tree*100
            ax2.plot(time_tree,errtree_spitzer,color='black')
        if show_hosoinu == True and plot_parts == False:
            ax1.plot(time_tree,hosoinurad_tree,color='grey',label='Hosokawa-Inutsuka solution')
            errtree_hosoinu = (ionrad_tree-hosoinurad_tree)/hosoinurad_tree*100
            ax2.plot(time_tree,errtree_hosoinu,color='grey')
        ax1.plot(time_tree,ionrad_tree,color='red',linewidth=1,label='SPH + MCRT with pseudo-particles')

    if plot_parts == True:
        
        file_exist = os.path.isfile('./time_rad_spitz_parts.txt')
        if file_exist == True:
            # read existing file
            time_parts,ionrad_parts,spitzrad_parts,hosoinurad_parts = np.loadtxt('time_rad_spitz_parts.txt',unpack=True)
        else:
            # compute ion-rad
            time_parts,ionrad_parts,spitzrad_parts,hosoinurad_parts = extract_data(False)
            # save to file 
            data = np.column_stack([time_parts,ionrad_parts,spitzrad_parts,hosoinurad_parts])
            np.savetxt('time_rad_spitz_parts.txt',data)
        
        ax1.plot(time_parts,ionrad_parts,color='blue',linewidth=1,label='SPH + MCRT with ind. particles')
        if show_spitzer == True:
            ax1.plot(time_parts,spitzrad_parts,color='black',label='Spitzer solution')
            errparts_spitzer = (ionrad_parts-spitzrad_parts)/spitzrad_parts*100
            ax2.plot(time_parts,errparts_spitzer,color='black')
        if show_hosoinu == True: 
            ax1.plot(time_parts,hosoinurad_parts,color='black',label='Hosokawa-Inutsuka solution')
            errparts_hosoinu = (ionrad_parts-hosoinurad_parts)/hosoinurad_parts*100
            ax2.plot(time_parts,errparts_hosoinu,color='grey')

    plt.rcParams.update({'font.size':14})
    ax1.legend(prop={'size': 11})
    ax1.set_xlabel('Time [Myr]',fontsize=14)
    ax1.set_ylabel('Ionization front radius [pc]',fontsize=14)
    ax1.tick_params(axis='x',labelsize=12)
    ax1.tick_params(axis='y',labelsize=12)
    ax1.minorticks_on()
    ax1.set_xlim([-0.005,0.125])
    
    ax2.set_xlabel('Time [Myr]',fontsize=14)
    ax2.set_ylabel('Residuals [%]',fontsize=11)
    ax2.set_xlim([-0.005,0.125])
    ax2.set_ylim([-9,4])
#    ax2.set_ylim([-6,17])
    ax2.hlines(y=0, xmin=-0.005, xmax=0.125, color='grey', linestyle='--', lw=1, alpha=0.3)
    ax2.tick_params(axis='x',labelsize=12)
    ax2.tick_params(axis='y',labelsize=12)
    ax2.minorticks_on()
    #plt.show()
    plt.savefig('ionization_front_radius_evol.svg')
    
    
def extract_data(tree):
    '''
    Returns time [Myr] vs ionization radius [pc] & spitzer radius [pc]
    '''
    time = []
    ionrad = []
    spitzrad = []
    hosoinurad = []
    
    pwd = os.getcwd()
    pardir = os.path.abspath(os.path.join(pwd, os.pardir))
    if tree == True:
        path = 'tree128_snapshots/nixyzhmf_*.txt'
    else:
        path = 'parts_snapshots/xyzhmf_*.txt'
    for filename in glob.glob(path):  
        print(filename)
        
        file = open(filename)
        for iline,line in enumerate(file):
            if (iline == 1):
                time_cgs = line 
                break
        time_cgs = float(time_cgs)

        nixyzhmf = np.loadtxt(filename,skiprows=3)
        if tree == True:
            x = nixyzhmf[:,2]
            y = nixyzhmf[:,3]
            z = nixyzhmf[:,4]
            nH = nixyzhmf[:,7]
        else:
            x = nixyzhmf[:,0]
            y = nixyzhmf[:,1]
            z = nixyzhmf[:,2]
            nH = nixyzhmf[:,5]
            
        rad = get_ionrad(x,y,z,nH)
        rad_cgs = rad * udist  
        
        # Calculate spitzer solution at time_cgs
        rad_strm = ((3.*lumin*mass_p**2)/(4.*np.pi*alpha*rhozero**2))**(1./3.)
        rad_spitz = rad_strm * (1. + 7./4.*(cs_ion*time_cgs/rad_strm))**(4./7.)
        
        # Calculate hosogawa-inutsuka solution 
        rad_hosoinu = rad_strm * (1. + 7./4.*np.sqrt(4./3.)*(cs_ion*time_cgs/rad_strm))**(4./7.)
        
        # Convert to Myr,pc
        time_myr = time_cgs/(1E6*365*24*60*60) 
        rad_pc = rad_cgs/3.08567758128E+18
        rad_spitz_pc = rad_spitz/3.08567758128E+18
        rad_hosoinu_pc = rad_hosoinu/3.08567758128E+18
        
        time.append(time_myr)
        ionrad.append(rad_pc)
        spitzrad.append(rad_spitz_pc)
        hosoinurad.append(rad_hosoinu_pc)
    
    ionrad = sort_array(time,ionrad)
    spitzrad = sort_array(time,spitzrad)
    hosoinurad = sort_array(time,hosoinurad)
    time = sort_array(time,time)

    return np.array(time), np.array(ionrad), np.array(spitzrad), np.array(hosoinurad)
   

def sort_array(array1,array2):

    array2_sorted = [x for _,x in sorted(zip(array1,array2))]

    return array2_sorted
    

def get_ionrad(x_points,y_points,z_points,nH_points):
    '''
    Find ionization front from one output file 
    Input: 
        - xyzf, node properties 
    Output:
        - ionrad, radius of ionization front [pc]
    '''
    
    r_points = np.sqrt(x_points**2+y_points**2+z_points**2)
    rmax = np.max(r_points) 
    
    # Sort onto 1-D grid of nH vs r 
    ncell = 5000
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
        
