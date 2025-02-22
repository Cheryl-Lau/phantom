
'''
Outflows through a shell around feedback source 
Option 1: plot the evolution of total energy/momentum through spherical surface  
Option 2: plot distribution of energy/momentum on spherical surface for a particular dumpfile   
'''

import numpy as np 
import matplotlib.pyplot as plt 
import glob 
import scipy.interpolate


def plot_tot_evol():

    time_evol = []
    kin_evol = []
    int_evol = []

    for filename in glob.glob('outflow_withhii_chnl02_*.dat'):  # cgs units 
        file = open(filename)

        kin_tot = 0
        int_tot = 0
        for iline,line in enumerate(file):
            line = line.split()
            if (iline == 1):
                time = line[0]
                time_evol.append(float(time)/(1E6*365*24*60*60))
            elif (iline >= 3): 
                x = line[0]
                y = line[1]
                rho = line[2]
                vel = line[3]
                u = line[4]
                thermpr = line[5]
                rampr = line[6]
                kin_tot += 0.5*float(rho)*float(vel)**2   # KE per unit volume 
                int_tot += float(rho)*float(u)            # IE per unit volume 

        kin_evol.append(kin_tot)
        int_evol.append(int_tot)

    plt.plot(time_evol,kin_evol)

    # plot cummulative 



def plot_surface():
    
    file = open('outflow_withhii_chnl02_01000.dat')

    x_map = []
    y_map = []
    rho_distri = []
    momen_distri = []
    kin_distri = []
    int_distri = []
    for iline,line in enumerate(file):
        line = line.split()
        if (iline == 1):
            time = line[0]
            time_Myr = float(time)/(1E6*365*24*60*60)
        elif (iline >= 3): 
            x = line[0]
            y = line[1]
            rho = line[2]
            vel = line[3]
            u = line[4]
            thermpr = line[5]
            rampr = line[6]
            x_map.append(float(x))
            y_map.append(float(y))
            rho_distri.append(float(rho))
            momen_distri.append(float(rho)*float(vel))       # momentum per unit volume 
            kin_distri.append(0.5*float(rho)*float(vel)**2)   # KE per unit volume 
            int_distri.append(float(rho)*float(u))            # IE per unit volume 

    longmin = -np.pi
    longmax = np.pi
    latmin = -np.pi/2 + 0.2
    latmax = np.pi/2 - 0.2

    fig1 = plt.figure(figsize=[7,3])
    ax1 = fig1.add_subplot(111)
    ax1.text(80,5,str(round(time_Myr,2))+' Myr')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    grid1 = np.log10(rho_distri).reshape(101,101)
    im1 = ax1.imshow(grid1,cmap='viridis',interpolation='spline16',extent=[longmin,longmax,latmin,latmax])
    fig1.colorbar(im1, orientation='vertical',label='density')

    fig2 = plt.figure(figsize=[7,3])
    ax2 = fig2.add_subplot(111)
    ax2.text(80,5,str(round(time_Myr,2))+' Myr')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    grid2 = np.log10(momen_distri).reshape(101,101)
    im2 = ax2.imshow(grid2,cmap='viridis',interpolation='spline16',extent=[longmin,longmax,latmin,latmax])
    fig2.colorbar(im2, orientation='vertical',label='momentum')

    fig3 = plt.figure(figsize=[7,3])
    ax3 = fig3.add_subplot(111)
    ax3.text(80,5,str(round(time_Myr,2))+' Myr')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    grid3 = np.log10(kin_distri).reshape(101,101)
    im3 = ax3.imshow(grid3,cmap='viridis',interpolation='spline16',extent=[longmin,longmax,latmin,latmax])
    fig3.colorbar(im3, orientation='vertical',label='kinetic energy')

    fig4 = plt.figure(figsize=[7,3])
    ax4 = fig4.add_subplot(111)
    ax4.text(80,5,str(round(time_Myr,2))+' Myr')
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    grid4 = np.log10(int_distri).reshape(101,101)
    im4 = ax4.imshow(grid4,cmap='viridis',interpolation='spline16',extent=[longmin,longmax,latmin,latmax])
    fig4.colorbar(im4, orientation='vertical',label='thermal energy')



def main(): 

#    plot_tot_evol()
    plot_surface()

    plt.show()
    

if __name__ == "__main__":
    main()
        









