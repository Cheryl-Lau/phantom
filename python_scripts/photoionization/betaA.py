#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 14:26:35 2023

@author: okitacscl1
"""

import numpy as np
import matplotlib.pyplot as plt 
import scipy as sp

x = [2500.,5000.,10000.,20000.]
y = [8.93E-13,5.42E-13,3.23E-13,1.88E-13]
x = np.array(x)
y = np.array(y)


def func(x,a,b):
#    return a* -x**(-1./2.) * np.exp(-b/x)
    return a* x**(-1./2.) + b



popt, pcov = sp.optimize.curve_fit(func,x,y)

print(*popt)

xsmooth = np.linspace(500,20000,100)
plt.plot(xsmooth, func(xsmooth,*popt))
plt.plot(x,y,'.')

plt.show()