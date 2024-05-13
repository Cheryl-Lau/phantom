# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt 
import numpy as np 
import scipy.optimize as sp

data = np.loadtxt('ms_lifetimes.txt',unpack=True) 

mass = data[0] 
time = data[1] * 1E6

plt.plot(mass,time,'o')
plt.show()

def fitfunc(mass,k,a):
    time = k*mass**(-a)
    return time


fitparam = sp.curve_fit(fitfunc,mass,time,p0=[1E6,2])

print(fitparam)

plt.plot(mass,fitfunc(mass,fitparam[0][0],fitparam[0][1]))


# For 8 M_sun 
time8 = fitfunc(8,fitparam[0][0],fitparam[0][1])
print('8Msun lifetime',time8)
 


































