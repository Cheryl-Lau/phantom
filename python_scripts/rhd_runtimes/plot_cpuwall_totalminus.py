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
    scale_cpu = grad_loc/grad_hpc

    grad_hpc, intercept_hpc = np.polyfit(sim_hpc, wall_hpc,1)
    grad_loc, intercept_loc = np.polyfit(sim_loc, wall_loc,1)
    scale_wall = grad_loc/grad_hpc

    return scale_cpu, scale_wall 


def correct_kennedy_speed_cmi():
    istep, cpu_hpc, wall_hpc = np.loadtxt('ten7particles_kennedy/set16/CMI_cpu_wall_time_record.txt', unpack=True)
    istep, cpu_loc, wall_loc = np.loadtxt('ten7particles/set5/CMI_cpu_wall_time_record.txt', unpack=True)

    cpumean_loc = np.mean(cpu_loc)
    cpumean_hpc = np.mean(cpu_hpc)
    scale_cmi_cpu = cpumean_loc/cpumean_hpc

    wallmean_loc = np.mean(wall_loc)
    wallmean_hpc = np.mean(wall_hpc)
    scale_cmi_wall = wallmean_loc/wallmean_hpc

    return scale_cmi_cpu, scale_cmi_wall


partnumsets = ['ten5particles','ten6particles','ten7particles','ten7particles_kennedy']
totpart = [144666,1426800,9938375,9938375]
colours = ['red','blue','black','grey']
labels = [r'$1.47\times10^5$'+' particles', r'$1.43\times10^6$'+' particles', \
            r'$0.99\times10^7$'+' particles', r'$0.99\times10^7$'+' particles on HPC']


dt = 6.701E-06

i = 0 


fig1, (ax2,ax1) = plt.subplots(nrows=2, sharex=True, subplot_kw=dict(frameon=True), figsize=(7,10), dpi=300)
plt.subplots_adjust(hspace=.0)

fig2, (ax4,ax3) = plt.subplots(nrows=2, sharex=True, subplot_kw=dict(frameon=True), figsize=(7,10), dpi=300)
plt.subplots_adjust(hspace=.0)


for numpartdir, allpart in zip(partnumsets, totpart):

    nppart_set = []
    cputime_set = []
    cpuerr_set = []
    cmitime_set = []
    walltime_set = []
    wallerr_set = []
    treetime_cpu_set = []
    treeerr_cpu_set = []
    treetime_wall_set = []
    treeerr_wall_set = []

    scale_fac_cpu, scale_fac_wall = correct_kennedy_speed()
    scale_fac_cmi_cpu, scale_fac_cmi_wall = correct_kennedy_speed_cmi()


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
            cpu_time = cpu_time * scale_fac_cpu
            wall_time = wall_time * scale_fac_wall


        # total cpu_time per step and error 
        grad_guess, intercept_guess = np.polyfit(sim_time, cpu_time,1)
        popt, pcov = curve_fit(linear_func, sim_time, cpu_time, np.array([grad_guess, intercept_guess]))
        grad_cpu = popt[0]
        graderr_cpu = np.sqrt(np.diag(pcov))[0]

        grad_guess, intercept_guess = np.polyfit(sim_time, wall_time,1)
        popt, pcov = curve_fit(linear_func, sim_time, wall_time, np.array([grad_guess, intercept_guess]))
        grad_wall = popt[0]
        graderr_wall = np.sqrt(np.diag(pcov))[0]


        # CMI wall time measured from phantom 
        cmifile_phantom = setname+"/CMI_cpu_wall_time_record.txt"
        istep, cpu_cmi, wall_cmi = np.loadtxt(cmifile_phantom, unpack=True)
        cpumean_cmi = np.mean(cpu_cmi)
        cpuerr_cmi = np.std(cpu_cmi)
        wallmean_cmi = np.median(wall_cmi)
        wallerr_cmi = np.std(wall_cmi)/np.sqrt(len(wall_cmi))

        # correct for difference in cpu performance 
        if (numpartdir == 'ten7particles_kennedy' or setname == 'ten7particles/set16_fromkennedy'):
            cpumean_cmi = cpumean_cmi * scale_fac_cmi_cpu
            cpuerr_cmi = cpuerr_cmi * scale_fac_cmi_cpu
            wallmean_cmi = wallmean_cmi * scale_fac_cmi_wall
            wallerr_cmi = wallerr_cmi * scale_fac_cmi_wall


        # Time difference
        if (numpartdir == 'ten5particles'):

            exclcmifile_phantom = setname+"/exclude_CMI_cpu_wall_time_record.txt"
            istep, cpu2, wall2, cpu1, wall1 = np.loadtxt(exclcmifile_phantom, unpack=True)

            cpudiff = []
            for istep in range(len(cpu2)-1):
                cpudiff.append(cpu1[istep+1] - cpu2[istep]) 
            cpudiff_avg = np.mean(cpudiff)
            cpudiff_err = np.std(cpudiff)

            walldiff = []
            for istep in range(len(wall2)-1):
                walldiff.append(wall1[istep+1] - wall2[istep]) 
            walldiff_avg = np.mean(walldiff)
            walldiff_err = np.std(walldiff)

        else: 
            cpudiff_avg = grad_cpu*dt - cpumean_cmi
            cpudiff_err = np.sqrt((graderr_cpu*dt)**2 + cpuerr_cmi**2)

            walldiff_avg = grad_wall*dt - wallmean_cmi
            walldiff_err = np.sqrt((graderr_wall*dt)**2 + wallerr_cmi**2)


        nppart_set.append(nppart)
        cputime_set.append(grad_cpu*dt)
        cpuerr_set.append(graderr_cpu*dt)
        cmitime_set.append(cpumean_cmi)
        walltime_set.append(grad_wall*dt)
        wallerr_set.append(graderr_wall*dt)
        treetime_cpu_set.append(cpudiff_avg)
        treeerr_cpu_set.append(cpudiff_err)
        treetime_wall_set.append(walldiff_avg)
        treeerr_wall_set.append(walldiff_err)


    # sort by nppart 
    nppart_sorted = np.sort(nppart_set)
    cputime_sorted = [i for _, i in sorted(zip(nppart_set, cputime_set))]
    cpuerr_sorted = [i for _, i in sorted(zip(nppart_set, cpuerr_set))]
    walltime_sorted = [i for _, i in sorted(zip(nppart_set, walltime_set))]
    wallerr_sorted = [i for _, i in sorted(zip(nppart_set, wallerr_set))]
    cmitime_sorted = [i for _, i in sorted(zip(nppart_set, cmitime_set))]
    treetime_cpu_sorted = [i for _, i in sorted(zip(nppart_set, treetime_cpu_set))]
    treeerr_cpu_sorted = [i for _, i in sorted(zip(nppart_set, treeerr_cpu_set))]
    treetime_wall_sorted = [i for _, i in sorted(zip(nppart_set, treetime_wall_set))]
    treeerr_wall_sorted = [i for _, i in sorted(zip(nppart_set, treeerr_wall_set))]


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


    # plot walltime vs pseudo-particles 
    ax2.errorbar(nppart_sorted[:-1], walltime_sorted[:-1], yerr=wallerr_sorted[:-1], fmt='s', markersize=2, elinewidth=1, color=colours[i], label=labels[i])

    # plot all-particles case 
    if (numpartdir != 'ten7particles'):
        ax2.errorbar(nppart_sorted[-1], walltime_sorted[-1], yerr=wallerr_sorted[-1], fmt='D', markersize=3, elinewidth=1, color=colours[i])

    # plot fit line
    nppart_curve = np.linspace(nppart_sorted.min(), nppart_sorted.max(), 300) 
    smoothed_curve = make_interp_spline(nppart_sorted, walltime_sorted, k=1) 
    walltime_curve = smoothed_curve(nppart_curve) 
    ax2.plot(nppart_curve, walltime_curve, linewidth=1, c=colours[i])
    


    # plot cputime on SPH side 
    ax3.errorbar(nppart_sorted[:-2], treetime_cpu_sorted[:-2], yerr=treeerr_cpu_sorted[:-2], fmt='s', markersize=2, elinewidth=1, color=colours[i], label=labels[i])

    # plot fit line
    smoothed_curve = make_interp_spline(nppart_sorted[:-2], treetime_cpu_sorted[:-2], k=1) 
    treetime_cpu_curve = smoothed_curve(nppart_curve[:-2]) 
    ax3.plot(nppart_curve[:-2], treetime_cpu_curve, linewidth=1, c=colours[i])


    # plot walltime on SPH side 
    ax4.errorbar(nppart_sorted[:-2], treetime_wall_sorted[:-2], yerr=treeerr_wall_sorted[:-2], fmt='s', markersize=2, elinewidth=1, color=colours[i], label=labels[i])

    # plot fit line
    smoothed_curve = make_interp_spline(nppart_sorted[:-1], treetime_wall_sorted[:-1], k=1) 
    treetime_wall_curve = smoothed_curve(nppart_curve[:-1]) 
    ax4.plot(nppart_curve[:-1], treetime_wall_curve, linewidth=1, c=colours[i])


    i += 1 


ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel('Total CPU time per step [s]', fontsize = 14)
ax1.set_xlabel('Number of pseudo-particles', fontsize = 14)
ax1.set_xlim([1.5e4,1.5e7])
ax1.set_ylim([7E0,8E4])
ax1.tick_params(axis='x', labelsize=12)
ax1.tick_params(axis='y', labelsize=12)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylabel('Total wall time per step [s]', fontsize = 14)
ax2.set_xlabel('Number of pseudo-particles', fontsize = 14)
ax2.set_xlim([1.5e4,1.5e7])
ax2.set_ylim([7E0,8E4])
ax2.legend()
ax2.tick_params(axis='x', labelsize=12)
ax2.tick_params(axis='y', labelsize=12)

fig1.savefig('total_time.png')


plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_ylabel('SPH-side CPU time per step [s]', fontsize = 14)
ax3.set_xlabel('Number of pseudo-particles', fontsize = 14)
ax3.set_xlim([1.5e4,1.5e7])
ax3.set_ylim([9e-1,9e4])
ax3.tick_params(axis='x', labelsize=12)
ax3.tick_params(axis='y', labelsize=12)

plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set_ylabel('SPH-side wall time per step [s]', fontsize = 14)
ax4.set_xlabel('Number of pseudo-particles', fontsize = 14)
ax4.set_xlim([1.5e4,1.5e7])
ax4.set_ylim([9e-1,9e4])
ax4.legend()
ax4.tick_params(axis='x', labelsize=12)
ax4.tick_params(axis='y', labelsize=12)

fig2.savefig('sph_side_time.png')











