import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import glob 
from scipy.optimize import curve_fit


partnumsets = ['ten5particles'] #,'ten6particles','ten7particles']
totpart = [97336] #,999999,9938375]
colours = ['red'] #,'blue','black']
labels = [r'$0.97\times10^5$'+' particles'] #, r'$1.43\times10^6$'+' particles',r'$0.99\times10^7$'+' particles']


def linear_func(x,k,c):
    return k*x + c


def get_avg_cputime_for_step(df,process):

    beforeprocess = 'before'+process
    afterprocess = 'after'+process 

    try:
        df_before = df[df['description'].str.contains(beforeprocess)]
        df_after = df[df['description'].str.contains(afterprocess)]

        cputime = df_after['cpu_elapsed'].to_numpy() - df_before['cpu_elapsed'].to_numpy()

    except: # for allparts runs 
        cputime = 0

    return np.mean(cputime), np.std(cputime)


def sort_array(array_ref,array_in):

    array_sorted = [i for _, i in sorted(zip(array_ref, array_in))]

    return array_sorted


def main(): 

    fig = plt.figure(figsize=[7,5])
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

            # number of pseudo-particles 
            nppartfile = setname+"/nppart"
            nppart = np.loadtxt(nppartfile, unpack=True)
            # number of cells 
            ncellfile = setname+"/ncell"
            ncell = np.loadtxt(ncellfile, unpack=True)

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
            df = pd.read_fwf(indruntimefile, sep=" ", header=None)
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


        ax1.errorbar(nppart_set,meancpu_set,yerr=stdcpu_set,marker='s',markersize=1,elinewidth=1,color=colours[c], label=labels[c])
        ax2.errorbar(nppart_set,meancpu_tree_set,yerr=stdcpu_tree_set,marker='s',markersize=1,elinewidth=1,color=colours[c], label=labels[c])
        ax3.errorbar(nppart_set,meancpu_hsol_set,yerr=stdcpu_hsol_set,marker='s',markersize=1,elinewidth=1,color=colours[c], label=labels[c])
        ax4.errorbar(nppart_set,meancpu_cmi_set,yerr=stdcpu_cmi_set,marker='s',markersize=1,elinewidth=1,color=colours[c], label=labels[c])

        c += 1

    ax1.set_ylabel('Total CPU time')
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax2.set_ylabel('Tree-walk CPU time')
    ax2.set_xscale('log')
    ax2.set_yscale('log')

    ax3.set_ylabel('h_solve CPU time')
    ax3.set_xscale('log')
    ax3.set_yscale('log')

    ax4.set_ylabel('CMI CPU time')
    ax4.set_xscale('log')
    ax4.set_yscale('log')


    plt.show()



main()




















