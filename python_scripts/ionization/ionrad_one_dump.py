# -*- coding: utf-8 -*-
"""
Plot ionization fronts for one snapshot 
"""
import numpy as np
import matplotlib.pyplot as plt 
import os
import math

# Params in cgs units 
ncell = 5000
alpha = 2.7E-13
lumin = 1E49
rhozero = 5.21E-21
mass_p = 1.67262158E-24
cs_ion = 12.85E5
udist = 3.086E+18

Myr = 1E6*365*24*60*60
pc = 3.08567758128E+18


def main():
    
    pwd = os.getcwd()
    pardir = os.path.abspath(os.path.join(pwd, os.pardir))
    filename1 = pardir+'\\tree128_snapshots\\nixyzhmf_00030.txt'
    filename2 = pardir+'\\tree128_snapshots\\nixyzhmf_00105.txt'
    
    fig = plt.figure(figsize=(8, 8),dpi=200)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    plt.subplots_adjust(hspace=.0)
    
    plot_snapshot(filename1,ax1)
    plot_snapshot(filename2,ax2)
    
    plt.savefig('ionfront_r0208.png')
    plt.show()
    
    return 


def plot_snapshot(filename,ax):
    
    time_Myr,x,y,z,nH,rad_spitz,rad_hosoinu = extract_data(filename)

    r = calculate_radius(x,y,z)
    r02,r08 = get_ionrad_0208(x,y,z,r,nH)
    
    ax.axvline(x=rad_spitz,linewidth=1,linestyle='-',color='black',label='Spitzer solution')
    ax.axvline(x=rad_hosoinu,linewidth=1,linestyle='-',color='grey',label='Hosokawa-Inutsuka solution')
    ax.axvline(x=r02,linewidth=1,linestyle='--',color='orange',label='$f_\mathrm{neu} = 0.2$')
    ax.axvline(x=r08,linewidth=1,linestyle='--',color='red',label='$f_\mathrm{neu} = 0.8$')

    ax.scatter(r,nH,s=0.1,color='royalblue',label='pseudo-particles')
    
    ax.text(0.05,0.95,'t = '+str(round(time_Myr,2))+' Myr',fontsize=12)
    ax.set_xlabel('Radius [pc]')
    ax.set_ylabel('Neutral fraction')
    ax.legend(fontsize=9,loc='lower right')
    ax.set_ylim([-0.05,1.1])
    ax.minorticks_on()
    
    return 


def extract_data(filename):

    time_cgs = np.loadtxt(filename,skiprows=1,max_rows=1)
    n,i,x,y,z,h,m,nH = np.loadtxt(filename,skiprows=3,unpack=True)

    # Calculate spitzer solution & H-I solution at time_cgs
    rad_strm = ((3.*lumin*mass_p**2)/(4.*np.pi*alpha*rhozero**2))**(1./3.)
    rad_spitz = rad_strm * (1. + 7./4.*(cs_ion*time_cgs/rad_strm))**(4./7.)
    rad_hosoinu = rad_strm * (1. + 7./4.*np.sqrt(4./3.)*(cs_ion*time_cgs/rad_strm))**(4./7.)
    
    # convert units 
    time_Myr = time_cgs/Myr
    rad_spitz_pc = rad_spitz/pc
    rad_hosoinu_pc = rad_hosoinu/pc
    
    return time_Myr,x,y,z,nH,rad_spitz_pc,rad_hosoinu_pc


def calculate_radius(x,y,z):
    r = np.sqrt(x**2 + y**2 + z**2)
    return r 
    

def get_ionrad_0208(x_points,y_points,z_points,r_points,nH_points):
    
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
    
    return r02,r08
    
    
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
    
    
    
    
    
    
    
    