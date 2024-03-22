import numpy as np
import matplotlib.pyplot as plt 
import os 
import glob 

partnumsets = ['ten5particles','ten6particles']
totpart = [14666,1426800]

for numpartdir, allpart in zip(partnumsets,totpart):

    nppart_set = []
    cputime_set = []

    path = numpartdir+"/set*"
    for setname in glob.glob(path): 

        # number of pseudo-particles 
        nppartfile = setname+"/nppart"
        nppart = np.loadtxt(nppartfile,unpack=True)

        # raw runtime 
        runtimefile = setname+"/cpu_wall_time_record.txt"
        istep, sim_time, cpu_time, wall_time = np.loadtxt(runtimefile,unpack=True)

        # cpu_time per sim_time
#        grad, intercept = np.polyfit(sim_time, cpu_time, 1)
        n = len(istep)-1
        grad = (cpu_time[n]-cpu_time[1])/(sim_time[n]-sim_time[1])

        if (nppart < 15000):
            print(setname)

        nppart_set.append(nppart)
        cputime_set.append(grad)


    # sort by nppart 
    nppart_sorted = np.sort(nppart_set)
    cputime_sorted = [x for _, x in sorted(zip(nppart_set, cputime_set))]
    plt.scatter(nppart_sorted,cputime_sorted)

plt.xscale('log')
plt.yscale('log')
plt.ylabel('CPU time per unit time in simulation [s]')
plt.xlabel('number of pseudo-particles')
plt.show()

        




