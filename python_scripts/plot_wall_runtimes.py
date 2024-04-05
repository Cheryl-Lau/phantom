import numpy as np
import matplotlib.pyplot as plt 
import glob 
from scipy.optimize import curve_fit
from scipy.interpolate import make_interp_spline 


def linear_func(x,k,c):
    return k*x + c


def correct_kennedy_speed():
    istep, sim_hpc, cpu_hpc, wall_hpc = np.loadtxt('ten7particles_kennedy/set16/cpu_wall_time_record.txt', unpack=True)
    istep, sim_loc, cpu_loc, wall_loc = np.loadtxt('ten7particles/set5/cpu_wall_time_record.txt', unpack=True)

    grad_hpc, intercept_hpc = np.polyfit(sim_hpc, wall_hpc,1)
    grad_loc, intercept_loc = np.polyfit(sim_loc, wall_loc,1)

    scale = grad_loc/grad_hpc
    return scale 


def correct_kennedy_speed_cmi():
    istep, cpu_hpc, wall_hpc = np.loadtxt('ten7particles_kennedy/set16/CMI_cpu_wall_time_record.txt', unpack=True)
    istep, cpu_loc, wall_loc = np.loadtxt('ten7particles/set5/CMI_cpu_wall_time_record.txt', unpack=True)

    wallmean_loc = np.mean(wall_loc)
    wallmean_hpc = np.mean(wall_hpc)

    scale_cmi = wallmean_loc/wallmean_hpc
    return scale_cmi


partnumsets = ['ten5particles','ten6particles','ten7particles','ten7particles_kennedy']
totpart = [144666,1426800,9938375,9938375]
colours = ['red','blue','black','grey']
labels = [r'$1.47\times10^5$'+' particles', r'$1.43\times10^6$'+' particles', \
            r'$0.99\times10^7$'+' particles', r'$0.99\times10^7$'+' particles on HPC']

dt = 6.701E-06

i = 0 

fig1 = plt.figure(dpi=200)
ax1 = fig1.add_subplot(111)
fig2 = plt.figure(dpi=200)
ax2 = fig2.add_subplot(111)


for numpartdir, allpart in zip(partnumsets, totpart):

    nppart_set = []
    walltime_set = []
    wallerr_set = []
    treetime_set = []
    treeerr_set = []

    scale_fac = correct_kennedy_speed()
    scale_fac_cmi = correct_kennedy_speed_cmi()

    # Hydro CPU time 
    hydrotimefile = 'noinject/'+numpartdir+'/cpu_wall_time_record.txt'
    istep, sim_hydro, cpu_hydro, wall_hydro = np.loadtxt(hydrotimefile, unpack=True)
    grad, intercept = np.polyfit(sim_hydro, wall_hydro, 1)
    hydrotime = grad*dt

    path = numpartdir+"/set*"
    for setname in glob.glob(path): 

        # number of pseudo-particles 
        nppartfile = setname+"/nppart"
        nppart = np.loadtxt(nppartfile, unpack=True)

        # raw runtime 
        runtimefile = setname+"/cpu_wall_time_record.txt"
        istep, sim_time, cpu_time, wall_time = np.loadtxt(runtimefile, unpack=True)

        # correct for difference in cpu performance 
        if (numpartdir == 'ten7particles_kennedy'):
            wall_time = wall_time * scale_fac

        # total cpu_time per step and error 
        grad_guess, intercept_guess = np.polyfit(sim_time, wall_time,1)
        popt, pcov = curve_fit(linear_func, sim_time, wall_time, np.array([grad_guess, intercept_guess]))
        grad = popt[0]
        graderr = np.sqrt(np.diag(pcov))[0]


        # CMI wall time measured from phantom 
        cmifile_phantom = setname+"/CMI_cpu_wall_time_record.txt"
        istep, cpu_cmi, wall_cmi = np.loadtxt(cmifile_phantom, unpack=True)
        wallmean_cmi = np.median(wall_cmi)
        wallerr_cmi = np.std(wall_cmi)/np.sqrt(len(wall_cmi))

        # correct for difference in cpu performance 
        if (numpartdir == 'ten7particles_kennedy'):
            wallmean_cmi = wallmean_cmi * scale_fac_cmi
            wallerr_cmi = wallerr_cmi * scale_fac_cmi


        # Time difference
        if (numpartdir == 'ten5particles'):

            exclcmifile_phantom = setname+"/exclude_CMI_cpu_wall_time_record.txt"
            istep, cpu2, wall2, cpu1, wall1 = np.loadtxt(exclcmifile_phantom, unpack=True)
            walldiff = []
            for istep in range(len(wall2)-1):
                walldiff.append(wall1[istep+1] - wall2[istep]) 
            walldiff_avg = np.mean(walldiff)
            walldiff_err = np.std(walldiff)

        else: 
            walldiff_avg = grad*dt - wallmean_cmi
            walldiff_err = np.sqrt((graderr*dt)**2 + wallerr_cmi**2)


        nppart_set.append(nppart)
        walltime_set.append(grad*dt)
        wallerr_set.append(graderr*dt)
        treetime_set.append(walldiff_avg)
        treeerr_set.append(walldiff_err)


    # sort by nppart 
    nppart_sorted = np.sort(nppart_set)
    walltime_sorted = [i for _, i in sorted(zip(nppart_set, walltime_set))]
    wallerr_sorted = [i for _, i in sorted(zip(nppart_set, wallerr_set))]
    treetime_sorted = [i for _, i in sorted(zip(nppart_set, treetime_set))]
    treeerr_sorted = [i for _, i in sorted(zip(nppart_set, treeerr_set))]


    # plot cputime vs pseudo-particles 
    ax1.errorbar(nppart_sorted[:-1], walltime_sorted[:-1], yerr=wallerr_sorted[:-1], fmt='s', markersize=2, elinewidth=1, color=colours[i], label=labels[i])

    # plot all-particles case 
    if (numpartdir != 'ten7particles'):
        ax1.errorbar(nppart_sorted[-1], walltime_sorted[-1], yerr=wallerr_sorted[-1], fmt='D', markersize=3, elinewidth=1, color=colours[i])

    # plot fit line
    nppart_curve = np.linspace(nppart_sorted.min(), nppart_sorted.max(), 300) 
    smoothed_curve = make_interp_spline(nppart_sorted, walltime_sorted, k=1) 
    walltime_curve = smoothed_curve(nppart_curve) 
    ax1.plot(nppart_curve, walltime_curve, linewidth=1, c=colours[i])
    

    # plot time on SPH side 
    ax2.errorbar(nppart_sorted[:-2], treetime_sorted[:-2], yerr=treeerr_sorted[:-2], fmt='s', markersize=2, elinewidth=1, color=colours[i], label=labels[i])

    # plot fit line
    smoothed_curve = make_interp_spline(nppart_sorted[:-2], treetime_sorted[:-2], k=1) 
    treetime_curve = smoothed_curve(nppart_curve[:-2]) 
    ax2.plot(nppart_curve[:-2], treetime_curve, linewidth=1, c=colours[i])

    i += 1 




ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel('Wall time per step [s]')
ax1.set_xlabel('Number of pseudo-particles')
ax1.set_xlim([1.5e4,1.5e7])
ax1.set_ylim([7E0,8E4])
ax1.legend()
fig1.savefig('wall_time_pseudoparts.png')

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylabel('SPH-side Wall time per step [s]')
ax2.set_xlabel('Number of pseudo-particles')
ax2.set_xlim([1.5e4,2e7])
#ax2.set_ylim([1E0,8E4])
ax2.legend()
fig2.savefig('sphside_wall_time_pseudoparts.png')


plt.show()

        




