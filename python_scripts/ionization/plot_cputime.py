# -*- coding: utf-8 -*-
"""
Compares the cpu runtime of tree and parts
"""

import numpy as np
import matplotlib.pyplot as plt 
import os 
import glob

pwd = os.getcwd()
pardir = os.path.abspath(os.path.join(pwd, os.pardir))

file_nodes = pardir+'\\timing_record\\cpu_wall_time_record_tree.txt'
file_parts = pardir+'\\timing_record\\cpu_wall_time_record_allparts.txt'

data_nodes = np.loadtxt(file_nodes)
data_parts = np.loadtxt(file_parts)

irun_nodes = data_nodes[:,0]
cpu_time_nodes = data_nodes[:,1]

irun_parts = data_parts[:,0]
cpu_time_parts = data_parts[:,1]

#
# Extract corresponding sim time from node ionization files 
#
sim_time_nodes_ext = []
cpu_time_nodes_ext = []

for irun in range(0,int(irun_nodes[-1]),10):
    
    path = pardir+'\\tree_snapshots\\nixyzhmf_*.txt'
    for filename in glob.glob(path):  
        
        label = ""
        for c in filename:
            if c.isdigit():
                label = label + c
        label = int(label)

        if (label == irun/10):     
            file = open(filename)
            for iline,line in enumerate(file):
                if (iline == 1):
                    time_cgs = line 
                    break
            time_cgs = float(time_cgs)
            
            sim_time_nodes_ext.append(time_cgs)
            cpu_time_nodes_ext.append(cpu_time_nodes[irun])
            
#
# Extract corresponding sim time from parts ionization files 
#
sim_time_parts_ext = []
cpu_time_parts_ext = []

for irun in range(0,int(irun_parts[-1]),10):
    
    path = pardir+'\\parts_snapshots\\xyzhmf_*.txt'
    for filename in glob.glob(path):  
        
        label = ""
        for c in filename:
            if c.isdigit():
                label = label + c
        label = int(label)
        
        if (label == irun/10):        
            file = open(filename)
            for iline,line in enumerate(file):
                if (iline == 1):
                    time_cgs = line 
                    break
            time_cgs = float(time_cgs)
            
            sim_time_parts_ext.append(time_cgs)
            cpu_time_parts_ext.append(cpu_time_parts[irun])
            
slope_parts,intercept = np.polyfit(cpu_time_parts_ext[50:],sim_time_parts_ext[50:],1)
slope_nodes,intercept = np.polyfit(cpu_time_nodes_ext[50:],sim_time_nodes_ext[50:],1)
print('parts',slope_parts)
print('nodes',slope_nodes)
print(slope_nodes/slope_parts)

plt.figure(figsize=[15,10],dpi=200)
plt.rcParams.update({'font.size': 20})
plt.plot(cpu_time_parts_ext,sim_time_parts_ext,linewidth=4,label='all ind. particles')
plt.plot(cpu_time_nodes_ext,sim_time_nodes_ext,linewidth=4,label='pseudo-particles')
plt.ylabel('Simulation time [s]',fontsize=18)
plt.xlabel('CPU time [s]',fontsize=18)
plt.legend()
plt.show()

plt.savefig('cpu_time_compare.png')