'''
For each density case, plot ionization front radius evolution for 
instant, implicit and explicit case, with analytical solutions 
'''

import numpy as np 
import matplotlib.pyplot as plt 
import glob
import os
import math

rho_dirs = ['rho21','rho20','rho19','rho18']
rhozeros = [1e-21,1e-20,1e-19,1e-18]
therm_dirs = ['instant','implicit','explicit']
labels = ['instant','implicit','explicit','Spitzer solution']
colours = ['red','blue','green','black']
texts = [r'$\rho=10^{-21}$'+r'$\mathrm{\ g \ cm^{-3}}$', r'$\rho=10^{-20}$'+r'$\mathrm{\ g \ cm^{-3}}$', r'$\rho=10^{-19}$'+r'$\mathrm{\ g \ cm^{-3}}$', r'$\rho=10^{-18}$'+r'$\mathrm{\ g \ cm^{-3}}$']

# Params in cgs units 
ncell = 1000
alpha = 2.7E-13
lumin = 1E49
mass_p = 1.67262158E-24
cs_ion = 12.85E5
udist = 3.086E+18

plot_spitzer = True


def main():

    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(15, 4))
    plt.subplots_adjust(hspace=0.5)

    # loop through densities 
    t = 0
    for rho_dir,rhozero,ax in zip(rho_dirs,rhozeros,axs.ravel()):
        c = 0 

        # loop through instant, implcit, explicit 
        maxr = 0
        for therm_dir in therm_dirs:
            path = rho_dir+"/"+therm_dir
            print(path)
            time, ionrad, spitzrad, hosoinurad = extract_data(path,rhozero)

            ax.plot(time,ionrad,linewidth=1.0,color=colours[c],label=labels[c])
            c += 1 

        if (plot_spitzer == True):
            ax.plot(time,spitzrad,'--',linewidth=1.0,color=colours[c],label=labels[c])


        ax.set_title(texts[t])
        ax.legend(loc = 'upper left')
        t += 1

    plt.savefig('ionrad.png')
    plt.show()

    return 



def extract_data(pardir,rhozero):
    '''
    Returns time [Myr] vs ionization radius [pc] & spitzer/hosogawa-inutsuka radius [pc]
    '''
    time = []
    ionrad = []
    spitzrad = []
    hosoinurad = []
    
    path = pardir+'/nixyzhmf_0*.txt'

    for filename in glob.glob(path):  
        
        time_cgs = np.loadtxt(filename,max_rows=1,skiprows=1)

        nixyzhmf = np.loadtxt(filename,skiprows=3)
        x = nixyzhmf[:,2]
        y = nixyzhmf[:,3]
        z = nixyzhmf[:,4]
        nH = nixyzhmf[:,7]

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
    
    time_sorted = sort_array(time,time)
    ionrad_sorted = sort_array(time,ionrad)
    spitzrad_sorted = sort_array(time,spitzrad)
    hosoinurad_sorted = sort_array(time,hosoinurad)

    return time_sorted, ionrad_sorted, spitzrad_sorted, hosoinurad_sorted
    
    
    
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
    

def sort_array(array_ref,array_in):

    array_sorted = [i for _, i in sorted(zip(array_ref, array_in))]

    return array_sorted



main()


