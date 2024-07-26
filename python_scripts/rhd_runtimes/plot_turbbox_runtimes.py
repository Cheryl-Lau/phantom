import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import glob 
from scipy.optimize import curve_fit

import warnings
warnings.filterwarnings('ignore')


partnumsets = ['ten5particles','ten6particles','ten7particles']
totpart = [97336,1000000,9938375]
colours = ['red','blue','black']
labels = [r'$0.97\times10^5$'+' particles', r'$1.00\times10^6$'+' particles',r'$0.99\times10^7$'+' particles']

partnumsets.reverse()  # plot ten5particles curve on top
totpart.reverse()
colours.reverse()
labels.reverse()


def linear_func(x,k,c):
    return k*x + c


def get_avg_cputime_for_step(df,process):

    beforeprocess = 'before'+process
    afterprocess = 'after'+process 

    try:
        df_before = df[df['description'].str.contains(beforeprocess)]
        df_after = df[df['description'].str.contains(afterprocess)]

        cputime = df_after['cpu_elapsed'].to_numpy() - df_before['cpu_elapsed'].to_numpy()
#        if (process=='_treewalk'):
#            print(cputime)

    except: # for allparts runs 
        cputime = 0

    return np.mean(cputime), np.std(cputime)


def sort_array(array_ref,array_in):

    array_sorted = [i for _, i in sorted(zip(array_ref, array_in))]

    return array_sorted


def main(): 

    fig = plt.figure(figsize=[9,6],dpi=200)

    ax1 = fig.add_subplot(221)  # total time 
    ax2 = fig.add_subplot(222)  # tree time 
    ax3 = fig.add_subplot(223)  # hsolve time
    ax4 = fig.add_subplot(224)  # cmi time 

    c = 0 

    for numpartdir, allpart in zip(partnumsets, totpart):

        # init storages 
        nppart_set = []
        meancpu_set,meancpu_tree_set,meancpu_hsol_set,meancpu_cmi_set = ([] for i in range(4))
        stdcpu_set,stdcpu_tree_set,stdcpu_hsol_set,stdcpu_cmi_set = ([] for i in range(4))

        # loop through each set 
        path = numpartdir+"/set*"
        for setname in glob.glob(path): 
#            print(setname)

            # number of pseudo-particles 
            nppartfile = setname+"/nppart"
            nppart = np.loadtxt(nppartfile, unpack=True)

            # read runtime file 
            runtimefile = setname+"/cpu_wall_time_record.txt"
            istep, sim_time, cpu_time, wall_time = np.loadtxt(runtimefile, unpack=True)

            # cpu_time per step 
            popt, pcov = curve_fit(linear_func, istep, cpu_time)
            meancpu = popt[0]
            stdcpu = np.sqrt(np.diag(pcov))[0]


            # read runtime file for individual processes 
            indruntimefile = setname+"/cpu_wall_time_indv_process.txt"
            # put into dataframe 
            df = pd.read_fwf(indruntimefile, widths = [5,30,20,20], header=None)
            df.columns = ["irun", "description", "cpu_elapsed", "wall_elapsed"]

            # CPU time for tree-walk  
            meancpu_tree, stdcpu_tree = get_avg_cputime_for_step(df,'_treewalk')
            # CPU time for h-solve 
            meancpu_hsol, stdcpu_hsol = get_avg_cputime_for_step(df,'_hsolve')
            # CPU time for collecting full set of pseduo-particles 
            meancpu_coll, stdcpu_coll = get_avg_cputime_for_step(df,'_collectnodes')
            # CPU time for CMI 
            meancpu_cmi, stdcpu_cmi = get_avg_cputime_for_step(df,'CMI')

            nppart_set.append(nppart)
            meancpu_set.append(meancpu)
            meancpu_tree_set.append(meancpu_tree)
            meancpu_hsol_set.append(meancpu_hsol)
            meancpu_cmi_set.append(meancpu_cmi)
            stdcpu_set.append(stdcpu)
            stdcpu_tree_set.append(stdcpu_tree)
            stdcpu_hsol_set.append(stdcpu_hsol)
            stdcpu_cmi_set.append(stdcpu_cmi)


        meancpu_set = sort_array(nppart_set,meancpu_set)
        meancpu_tree_set = sort_array(nppart_set,meancpu_tree_set)
        meancpu_hsol_set = sort_array(nppart_set,meancpu_hsol_set)
        meancpu_cmi_set = sort_array(nppart_set,meancpu_cmi_set)
        stdcpu_tree_set = sort_array(nppart_set,stdcpu_tree_set)
        stdcpu_hsol_set = sort_array(nppart_set,stdcpu_hsol_set)
        stdcpu_cmi_set = sort_array(nppart_set,stdcpu_cmi_set)
        nppart_set = sort_array(nppart_set,nppart_set)

        ax1.errorbar(nppart_set,meancpu_set,yerr=stdcpu_set,fmt='s',markersize=2,elinewidth=1,color=colours[c],label=labels[c])
        ax1.scatter(nppart_set[-1],meancpu_set[-1],marker='*',s=20,color=colours[c])
        ax2.errorbar(nppart_set,meancpu_tree_set,yerr=stdcpu_tree_set,fmt='s',markersize=2,elinewidth=1,color=colours[c],label=labels[c])
        ax3.errorbar(nppart_set,meancpu_hsol_set,yerr=stdcpu_hsol_set,fmt='s',markersize=2,elinewidth=1,color=colours[c],label=labels[c])
        ax4.errorbar(nppart_set,meancpu_cmi_set,yerr=stdcpu_cmi_set,fmt='s',markersize=2,elinewidth=1,color=colours[c],label=labels[c])
        ax4.scatter(nppart_set[-1],meancpu_cmi_set[-1],marker='*',s=20,color=colours[c])

        c += 1

    ax1.set_ylabel('CPU time [s]')
    ax1.set_xlabel('Number of pseudo-particles')
    ax1.set_xlim([2.7e3,2e7])
    ax1.set_ylim([5e2,1.5e5])
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.text(5e3,6e4,'Total')
    ax1.legend(loc='lower right',fontsize=8)

    ax2.set_ylabel('CPU time [s]')
    ax2.set_xlabel('Number of pseudo-particles')
    ax2.set_xlim([2.7e3,2e7])
    ax2.set_ylim([8e-3,2e1])
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.text(5e3,6e0,'Tree-walk')
    ax2.legend(loc='lower right',fontsize=8)

    ax3.set_ylabel('CPU time [s]')
    ax3.set_xlabel('Number of pseudo-particles')
    ax3.set_xlim([2.7e3,2e7])
    ax3.set_ylim([2e-1,2e5])
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.text(5e3,3e4,'Smoothing length iterate')
    ax3.text(5e3,8e3,'(with Neighbour-find)')
    ax3.legend(loc='lower right',fontsize=8)

    ax4.set_ylabel('CPU time [s]')
    ax4.set_xlabel('Number of pseudo-particles')
    ax4.set_xlim([2.7e3,2e7])
    ax4.set_ylim([5e2,2e5])
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.text(5e3,8e4,'Density-mapping + MCRT')
    ax4.legend(loc='lower right',fontsize=8)


    fig.tight_layout()
    plt.savefig('turbbox_runtimes.png')
    plt.show()



main()




















