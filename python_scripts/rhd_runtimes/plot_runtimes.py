import numpy as np
import matplotlib.pyplot as plt 
import glob 
from scipy.optimize import curve_fit


def linear_func(x,k,c):
    return k*x + c


partnumsets = ['ten5particles','ten6particles','ten7particles','ten7particles_kennedy']
totpart = [14666,1426800,9938375,9938375]
colours = ['red','blue','black','grey']
labels = [r'$1.47\times10^5$'+' particles', r'$1.43\times10^6$'+' particles', \
            r'$0.99\times10^7$'+' particles', r'$0.99\times10^7$'+' particles on HPC']

c = 0 

for numpartdir, allpart in zip(partnumsets, totpart):

    nppart_set = []
    cputime_set = []
    graderr_set = []

    path = numpartdir+"/set*"
    for setname in glob.glob(path): 

        # number of pseudo-particles 
        nppartfile = setname+"/nppart"
        nppart = np.loadtxt(nppartfile, unpack=True)

        # raw runtime 
        runtimefile = setname+"/cpu_wall_time_record.txt"
        istep, sim_time, cpu_time, wall_time = np.loadtxt(runtimefile, unpack=True)

        # cpu_time per sim_time and error 
        popt, pcov = curve_fit(linear_func, sim_time, cpu_time)
        grad = popt[0]
        graderr = np.sqrt(np.diag(pcov))[0]

        nppart_set.append(nppart)
        cputime_set.append(grad)
        graderr_set.append(graderr)

    # sort by nppart 
    nppart_sorted = np.sort(nppart_set)
    cputime_sorted = [i for _, i in sorted(zip(nppart_set, cputime_set))]
    graderr_sorted = [i for _, i in sorted(zip(nppart_set, graderr_set))]

    plt.errorbar(nppart_sorted, cputime_sorted, yerr=graderr_sorted, fmt='.', color=colours[c], label=labels[c])
    c += 1 

plt.xscale('log')
plt.yscale('log')
plt.ylabel('CPU time per unit time in simulation [s]')
plt.xlabel('number of pseudo-particles')
plt.legend()
plt.show()

        




