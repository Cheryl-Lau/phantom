'''
Plotting the increase in ionized mass with different tree resol param K  
'''
import numpy as np 
import matplotlib.pyplot as plt 
import glob 

resol_dirs = ['K_30','K_100','K_500']
colours = ['red','green','blue']
labels = ['K = 30','K = 100','K = 500']

Myr = 1e6*365*24*60*60 

def main():
    
    fig = plt.figure(figsize=[6,4],dpi=200)
    ax1 = fig.add_subplot(111)

    # loop through K 
    for c,resol_dir in enumerate(resol_dirs):

        time_evol = []
        mass_evol = []

        path = resol_dir+"/ionized_mass_turbbox_*.dat"
        for dumpfile in glob.glob(path): 
            print(dumpfile)

            time = np.loadtxt(dumpfile,max_rows=1,skiprows=1)
            ionized_mass = np.loadtxt(dumpfile,max_rows=3,skiprows=3)

            time_evol.append(time/Myr)
            mass_evol.append(ionized_mass) 

        time_sorted = sort_array(time_evol,time_evol)
        mass_sorted = sort_array(time_evol,mass_evol)

        ax1.plot(time_sorted,mass_sorted,color=colours[c],label=labels[c])

    ax1.set_xlabel('Time [Myr]')
    ax1.set_ylabel('Ionized mass ['+r'$\mathrm{M}_{\odot}$'+']')
    ax1.legend()
    plt.savefig('ionized_mass_treeresol.png')
    plt.show()

    return 


def sort_array(array_ref,array_in):

    array_sorted = [i for _, i in sorted(zip(array_ref, array_in))]

    return array_sorted


main()



