#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

data = np.loadtxt('photon_flux_diazmiller.txt',skiprows=2)

solarm = 1.9891E33  

fit_mass = False
fit_temp = False

temp = data[:,0]
mass_solarm = data[:,1]
mass = mass_solarm * solarm
flux = data[:,2]   # THIS IS LOG10(FLUX)
plt.scatter(mass,temp)



def mass_func(x, a, c):
    return a*np.log10(x) + c


mass_sep = 30.0 * solarm 

xaxis = np.linspace(mass[0],mass[-1],100)

ihighmass = np.where(mass > mass_sep)[0]
popt, pcov = curve_fit(mass_func,mass[ihighmass],flux[ihighmass])
if (fit_mass):
    plt.plot(xaxis,mass_func(xaxis,*popt),'r-')
print('high mass reg',*popt)


ilowmass = np.where(mass < mass_sep)[0]
popt, pcov = curve_fit(mass_func,mass[ilowmass],flux[ilowmass])
if (fit_mass):
    plt.plot(xaxis,mass_func(xaxis,*popt),'g-')
print('low mass reg',*popt)


def temp_func(x,a,b,c,d):
    return a*np.exp(b*x + c) + d

xaxis = np.linspace(flux[0],flux[-1]-2,100)

popt, pcov = curve_fit(temp_func,flux,temp)
if (fit_temp): 
    plt.plot(xaxis,temp_func(xaxis,*popt))
print('q-temp',*popt)







plt.show()

