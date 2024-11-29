

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit


rad = [6e-2,8e-2,1e-1,2e-1,4e-1,6e-1,8e-1,1e0,2e0,4e0]
pc = 3.086E+18


def constant_func(r,c):
    return r**0 * c

def powerlaw_func(r,k,c):
    return r**k * c 


def plot_rho_sizescale(ax,filename,colourstr,labelstr): 

    df = pd.read_fwf(filename)
    df.columns = ['ip','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10']

    r_arr = []
    rho_avg = []
    rho_err = []
    for irad in range(len(rad)):
        r = rad[irad]
        colname = 'r'+str(irad+1)
        rhos = df[colname]
        r_arr.append(r)
        rho_avg.append(np.mean(rhos))
        rho_err.append(np.std(rhos))
        
    ax.errorbar(r_arr,rho_avg,yerr=rho_err,fmt='s',color=colourstr,markersize=5,label=labelstr)

    return r_arr,rho_avg,rho_err


def plot_colrho_sizescale(ax,filename,colourstr,labelstr): 

    df = pd.read_fwf(filename)
    df.columns = ['ip','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10']

    r_arr = []
    colrho_avg = []
    colrho_err = []
    for irad in range(len(rad)):
        r = rad[irad]
        colname = 'r'+str(irad+1)
        rhos = df[colname]
        r_arr.append(r)
        colrho_avg.append(np.mean(rhos)*(r*pc))
        colrho_err.append(np.std(rhos)*(r*pc))
        
    ax.errorbar(r_arr,colrho_avg,yerr=colrho_err,fmt='s',color=colourstr,markersize=5,label=labelstr)

    return r_arr,colrho_avg,colrho_err


def fit_curve(ax,func,r_arr,prop_vals,guess,colourstr,labelstr):

    r_val = np.logspace(np.log10(rad[0]*0.8),np.log10(rad[-1]*1.1),100)
    popt,pcov = curve_fit(func,r_arr,prop_vals,p0=guess)
    ax.plot(r_val,func(r_val,*popt),color=colourstr,label=labelstr)

    return popt



fig, (ax1,ax2) = plt.subplots(nrows=2, sharey=False, subplot_kw=dict(frameon=True), figsize=(5,7))


### Density

# Initial 
r_arr,rho_avg,rho_err = plot_rho_sizescale(ax1,'density_sizescale_cloud_20_10_00100.dat','grey','initial')
popt = fit_curve(ax1,constant_func,r_arr,rho_avg,(1e-11),'pink','initial '+r'$\rho = 10^{-10.4}\ R^{0}$')
print(*popt)

# After photoion
r_arr,rho_avg,rho_err = plot_rho_sizescale(ax1,'density_sizescale_cloud_20_10_clsink139_v2_06456.dat','black','after ion')
popt = fit_curve(ax1,powerlaw_func,r_arr,rho_avg,(-0.8,1e-20),'red','after ion '+r'$\rho = 10^{-19.8}\ R^{-0.8}$')
print(*popt)

# Only ionized gas 
r_arr,rho_avg,rho_err = plot_rho_sizescale(ax1,'density_sizescale_ionized_cloud_20_10_clsink139_v2_06456.dat','navy','ionized gas only')
popt = fit_curve(ax1,powerlaw_func,r_arr,rho_avg,(-3.0,1e-20),'royalblue','ionized '+r'$\rho = 10^{-22.0}\ R^{-2.0}$')
print(*popt)

ax1.set_xlabel('size-scale [pc]')
ax1.set_ylabel('density [$\mathrm{g\ cm^{-3}}$]')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylim([3e-40,2e-6])
ax1.legend(loc='lower right',fontsize=6)


### Column density

# Initial 
r_arr,colrho_avg,colrho_err = plot_colrho_sizescale(ax2,'density_sizescale_cloud_20_10_00100.dat','grey','initial')
popt = fit_curve(ax2,powerlaw_func,r_arr,colrho_avg,(1.0,2.46),'pink','initial '+r'$\rho R = 10^{8.08}\ R^{1.05}$')
print(*popt)

# After photoion
r_arr,colrho_avg,colrho_err = plot_colrho_sizescale(ax2,'density_sizescale_cloud_20_10_clsink139_v2_06456.dat','black','after ion')
popt = fit_curve(ax2,constant_func,r_arr,colrho_avg,(2.46),'red','after ion '+r'$\rho R = 0.05\ R^{0}$')
print(*popt)

# Only ionized gas 
r_arr,colrho_avg,colrho_err = plot_colrho_sizescale(ax2,'density_sizescale_ionized_cloud_20_10_clsink139_v2_06456.dat','navy','ionized gas only')
popt = fit_curve(ax2,powerlaw_func,r_arr,colrho_avg,(-2.0,2.46),'royalblue','ionized '+r'$\rho R = 10^{-3.68}\ R^{-1.20}$')
print(*popt)



ax2.set_xlabel('size-scale [pc]')
ax2.set_ylabel('column density [$\mathrm{g\ cm^{-2}}$]')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim([3e-17,1.8e10])
ax2.legend(loc='lower right',fontsize=6)



fig.tight_layout(pad=0.5)
plt.savefig('density_sizescale.png',dpi=200)

plt.show()




