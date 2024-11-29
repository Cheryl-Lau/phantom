

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit


rad = [6e-2,8e-2,1e-1,2e-1,4e-1,6e-1,8e-1,1e0,2e0,4e0]


def sigma_func2(r,c):
    return r**(1/2) * c 

def sigma_func3(r,c):
    return r**(1/3) * c 


def plot_sigma_relation(ax,filename,colourstr,labelstr):

    df = pd.read_fwf(filename)
    df.columns = ['ip','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10']

    r_arr = []
    sigma_avg = []
    sigma_err = []
    for irad in range(len(rad)):
        r = rad[irad]
        colname = 'r'+str(irad+1)
        sigmas = df[colname].values
        #ax.scatter([r]*len(sigmas),sigmas,s=1,c=colourstr,marker='.')
        r_arr.append(r)
        sigma_avg.append(np.mean(sigmas))
        sigma_err.append(np.std(sigmas))

    ax.errorbar(r_arr,sigma_avg,yerr=sigma_err,fmt='s',color=colourstr,markersize=5,label=labelstr)

    return r_arr,sigma_avg,sigma_err


def fit_curve(func,r_arr,sigma_avg,colourstr,labelstr):

    r_val = np.logspace(np.log10(rad[0]*0.8),np.log10(rad[-1]*1.1),100)
    popt,pcov = curve_fit(func,r_arr,sigma_avg)
    ax.plot(r_val,func(r_val,*popt),color=colourstr,label=labelstr)

    return popt


fig, ax = plt.subplots(ncols=1, sharey=False, subplot_kw=dict(frameon=True), figsize=(5.5,4))

# Initial 
r_arr,sigma_avg,sigma_err = plot_sigma_relation(ax,'velocity_dispersion_cloud_20_10_00100.dat','grey','initial')
popt = fit_curve(sigma_func2,r_arr,sigma_avg,'pink','initial $\sigma_v = 10^{5.17}\ R^{1/2}$')
print(np.log10(*popt))


# After photoion 
r_arr,sigma_avg,sigma_err = plot_sigma_relation(ax,'velocity_dispersion_cloud_20_10_clsink139_v2_06456.dat','black','after ion')
popt = fit_curve(sigma_func2,r_arr,sigma_avg,'red','after ion $\sigma_v = 10^{5.26}\ R^{1/2}$')
print(np.log10(*popt))


# Only ionized particles 
r_arr,sigma_avg,sigma_err = plot_sigma_relation(ax,'velocity_dispersion_ionized_cloud_20_10_clsink139_v2_06456.dat','navy','ionized gas only')
sigma_avg.pop(-1)
sigma_avg.pop(-1)
r_arr.pop(-1)
r_arr.pop(-1)
popt = fit_curve(sigma_func2,r_arr,sigma_avg,'royalblue','ionized $\sigma_v = 10^{4.73}\ R^{1/2}$')
print(np.log10(*popt))
popt = fit_curve(sigma_func3,r_arr,sigma_avg,'aqua','ionized $\sigma_v = 10^{4.70}\ R^{1/3}$')
print(np.log10(*popt))




ax.set_xlabel('size-scale [pc]')
ax.set_ylabel('velocity dispersion [$\mathrm{cm\ s^{-1}}$]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([3e2,1.8e6])
ax.legend(loc='lower right',fontsize=8)

fig.tight_layout(pad=0.5)
plt.savefig('vel_dispersion.png',dpi=200)

plt.show()

