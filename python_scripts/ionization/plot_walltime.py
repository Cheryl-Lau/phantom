# -*- coding: utf-8 -*-
"""
Compares the wall runtime of tree and parts
"""

import numpy as np
import matplotlib.pyplot as plt 
import os 

pwd = os.getcwd()
pardir = os.path.abspath(os.path.join(pwd, os.pardir))

file_nodes = pardir+'\\timing_record\\cpu_wall_time_record_tree2.txt'
file_parts = pardir+'\\timing_record\\cpu_wall_time_record_allparts.txt'

data_nodes = np.loadtxt(file_nodes)
data_parts = np.loadtxt(file_parts)

irun_nodes = data_nodes[:,0]
wall_time_nodes = data_nodes[:,2]

irun_parts = data_parts[:,0]
wall_time_parts = data_parts[:,2]

slope_parts,intercept = np.polyfit(wall_time_parts,irun_parts,1)
slope_nodes,intercept = np.polyfit(wall_time_nodes,irun_nodes,1)
print('parts',slope_parts)
print('nodes',slope_nodes)
print(slope_nodes/slope_parts)

plt.plot(wall_time_parts,irun_parts,label='around particles')
plt.plot(wall_time_nodes,irun_nodes,label='around nodes')
plt.xlabel('Wall time [s]')
plt.ylabel('Number of CMI-calls')
plt.legend()
plt.show()

plt.savefig('wall_time_compare.png')