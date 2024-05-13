#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 18:11:45 2023

@author: okitacscl1
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('gamma_nH_distri.txt')

nH = data[:,0]
gamma = data[:,1] 

plt.plot(nH,gamma)
plt.xlabel('nH')
plt.ylabel('gamma [erg/s]')
plt.yscale('log')
plt.show()