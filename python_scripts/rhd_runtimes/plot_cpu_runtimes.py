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

    grad_hpc, intercept_hpc = np.polyfit(sim_hpc, cpu_hpc,1)
    grad_loc, intercept_loc = np.polyfit(sim_loc, cpu_loc,1)

    scale = grad_loc/grad_hpc
    return scale 


def correct_kennedy_speed_cmi():
    istep, cpu_hpc, wall_hpc = np.loadtxt('ten7particles_kennedy/set16/CMI_cpu_wall_time_record.txt', unpack=True)
    istep, cpu_loc, wall_loc = np.loadtxt('ten7particles/set5/CMI_cpu_wall_time_record.txt', unpack=True)

    cpumean_loc = np.mean(cpu_loc)
    cpumean_hpc = np.mean(cpu_hpc)

    scale_cmi = cpumean_loc/cpumean_hpc
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
fig3 = plt.figure(dpi=200)
ax3 = fig3.add_subplot(111)


for numpartdir, allpart in zip(partnumsets, totpart):

    nppart_set = []
    cputime_set = []
    cpuerr_set = []
    cmitime_set = []
    treetime_set = []
    treeerr_set = []

    scale_fac = correct_kennedy_speed()
    scale_fac_cmi = correct_kennedy_speed_cmi()

    # Hydro CPU time 
    hydrotimefile = 'noinject/'+numpartdir+'/cpu_wall_time_record.txt'
    istep, sim_hydro, cpu_hydro, wall_hydro = np.loadtxt(hydrotimefile, unpack=True)
    grad, intercept = np.polyfit(sim_hydro, cpu_hydro, 1)
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
        if (numpartdir == 'ten7particles_kennedy' or setname == 'ten7particles/set16_fromkennedy'):
            cpu_time = cpu_time * scale_fac

        # total cpu_time per step and error 
        grad_guess, intercept_guess = np.polyfit(sim_time, cpu_time,1)
        popt, pcov = curve_fit(linear_func, sim_time, cpu_time, np.array([grad_guess, intercept_guess]))
        grad = popt[0]
        graderr = np.sqrt(np.diag(pcov))[0]


        # CMI cpu time measured from phantom 
        cmifile_phantom = setname+"/CMI_cpu_wall_time_record.txt"
        istep, cpu_cmi, wall_cmi = np.loadtxt(cmifile_phantom, unpack=True)
        cpumean_cmi = np.mean(cpu_cmi)
        cpuerr_cmi = np.std(cpu_cmi)

        # correct for difference in cpu performance 
        if (numpartdir == 'ten7particles_kennedy' or setname == 'ten7particles/set16_fromkennedy'):
            cpumean_cmi = cpumean_cmi * scale_fac_cmi
            cpuerr_cmi = cpuerr_cmi * scale_fac_cmi


        # Time difference
        if (numpartdir == 'ten5particles'):

            exclcmifile_phantom = setname+"/exclude_CMI_cpu_wall_time_record.txt"
            istep, cpu2, wall2, cpu1, wall1 = np.loadtxt(exclcmifile_phantom, unpack=True)
            cpudiff = []
            for istep in range(len(cpu2)-1):
                cpudiff.append(cpu1[istep+1] - cpu2[istep]) 
            cpudiff_avg = np.mean(cpudiff)
            cpudiff_err = np.std(cpudiff)

        else: 
            cpudiff_avg = grad*dt - cpumean_cmi
            cpudiff_err = np.sqrt((graderr*dt)**2 + cpuerr_cmi**2)

        nppart_set.append(nppart)
        cputime_set.append(grad*dt)
        cpuerr_set.append(graderr*dt)
        cmitime_set.append(cpumean_cmi)
        treetime_set.append(cpudiff_avg)
        treeerr_set.append(cpudiff_err)

    # sort by nppart 
    nppart_sorted = np.sort(nppart_set)
    cputime_sorted = [i for _, i in sorted(zip(nppart_set, cputime_set))]
    cpuerr_sorted = [i for _, i in sorted(zip(nppart_set, cpuerr_set))]
    cmitime_sorted = [i for _, i in sorted(zip(nppart_set, cmitime_set))]
    treetime_sorted = [i for _, i in sorted(zip(nppart_set, treetime_set))]
    treeerr_sorted = [i for _, i in sorted(zip(nppart_set, treeerr_set))]


    # plot cputime vs pseudo-particles 
    ax1.errorbar(nppart_sorted[:-1], cputime_sorted[:-1], yerr=cpuerr_sorted[:-1], fmt='s', markersize=2, elinewidth=1, color=colours[i], label=labels[i])

    # plot all-particles case 
    if (numpartdir != 'ten7particles'):
        ax1.errorbar(nppart_sorted[-1], cputime_sorted[-1], yerr=cpuerr_sorted[-1], fmt='D', markersize=3, elinewidth=1, color=colours[i])

    # plot fit line
    nppart_curve = np.linspace(nppart_sorted.min(), nppart_sorted.max(), 300) 
    smoothed_curve = make_interp_spline(nppart_sorted, cputime_sorted, k=1) 
    cputime_curve = smoothed_curve(nppart_curve) 
    ax1.plot(nppart_curve, cputime_curve, linewidth=1, c=colours[i])


    # plot CMI time
    ax2.scatter(nppart_sorted, cmitime_sorted, s=2, color=colours[i], label=labels[i])


    # plot time on SPH side 
    ax3.errorbar(nppart_sorted[:-2], treetime_sorted[:-2], yerr=treeerr_sorted[:-2], fmt='s', markersize=2, elinewidth=1, color=colours[i], label=labels[i])

    # plot fit line
    smoothed_curve = make_interp_spline(nppart_sorted[:-2], treetime_sorted[:-2], k=1) 
    treetime_curve = smoothed_curve(nppart_curve[:-1]) 
    ax3.plot(nppart_curve[:-1], treetime_curve, linewidth=1, c=colours[i])


    i += 1 


ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel('Total CPU time per step [s]')
ax1.set_xlabel('Number of pseudo-particles')
ax1.set_xlim([1.5e4,1.5e7])
ax1.set_ylim([7E0,8E4])
ax1.legend()
fig1.savefig('cpu_time_pseudoparts.png')

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylabel('CMI CPU time per step [s]')
ax2.set_xlabel('Number of pseudo-particles')
ax2.set_xlim([1.5e4,2e7])
ax2.legend()
fig2.savefig('cmi_cpu_time_pseudoparts.png')


ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_ylabel('SPH-side CPU time per step [s]')
ax3.set_xlabel('Number of pseudo-particles')
ax3.set_xlim([1.5e4,1.5e7])
ax3.set_ylim([9e-1,9e4])
ax3.legend()
fig3.savefig('sphside_cpu_time_pseudoparts.png')


plt.show()

        




