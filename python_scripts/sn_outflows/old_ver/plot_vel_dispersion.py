

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit


def sigma_func(r,c):
    return r**(1/2) * c 


df = pd.read_fwf('velocity_dispersion_cloud_20_10_clsink139_v2_06456.dat')

df.columns = ['ip','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10']
rad = [6e-2,8e-2,1e-1,2e-1,4e-1,6e-1,8e-1,1e0,2e0,4e0]



fig, ax = plt.subplots(ncols=1, sharey=False, subplot_kw=dict(frameon=True), figsize=(5.5,4))


r_arr = []
sigma_avg = []
sigma_err = []
for irad in range(len(rad)):
    r = rad[irad]
    colname = 'r'+str(irad+1)
    sigmas = df[colname]
    r_arr.append(r)
    sigma_avg.append(np.mean(sigmas))
    sigma_err.append(np.std(sigmas))
    
ax.errorbar(r_arr,sigma_avg,yerr=sigma_err,fmt='s',color='black',markersize=5)


# Fit line 
r_val = np.logspace(np.log10(rad[0]*0.8),np.log10(rad[-1]*1.1),100)
popt,pcov = curve_fit(sigma_func,r_arr,sigma_avg)
print('fitted param',*popt)
ax.plot(r_val,sigma_func(r_val,*popt),color='red',label='$\sigma_v = 10^{5.4}\ R^{1/2}$')


ax.set_xlabel('size-scale [pc]')
ax.set_ylabel('velocity dispersion [$\mathrm{cm\ s^{-1}}$]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([8e3,1.8e6])
ax.legend(loc='lower right')

fig.tight_layout(pad=0.5)
plt.savefig('vel_dispersion.png',dpi=200)

plt.show()

