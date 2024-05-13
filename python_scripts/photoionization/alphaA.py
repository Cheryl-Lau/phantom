#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:17:47 2023

@author: okitacscl1
"""

import numpy as np
import matplotlib.pyplot as plt 
from scipy import optimize

x = [5000.,10000.,20000.]
y = [6.82E-13,4.18E-13,2.51E-13]
x = np.array(x)
y = np.array(y)


def func(x,a,b):
#    return a* -x**(-1./2.) * np.exp(-b/x)
    return a* x**(-1./2.) + b



popt, pcov = optimize.curve_fit(func,x,y)

print(*popt)

xsmooth = np.linspace(500,20000,100)
plt.plot(xsmooth, func(xsmooth,*popt))
plt.plot(x,y,'.')

plt.show()