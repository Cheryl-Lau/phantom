
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit



rad = [6e-2,8e-2,1e-1,2e-1,4e-1,6e-1,8e-1,1e0,2e0,4e0]


def constant_func(r,c):
    return r**0 * c

def powerlaw_func(r,k,c):
    return r**k * c 


def plot_virial(ax,filename,colourstr,labelstr):

    df_vir = pd.read_fwf(filename)
    df_vir.columns = ['ip','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10']

    r_arr = []
    virial_avg = []
    virial_err = []
    for irad in range(len(rad)):
        r = rad[irad]
        print(r)
        colname = 'r'+str(irad+1)
        virials = df_vir[colname].astype(float)
        posvirials = virials.loc[virials > 0]
        r_arr.append(r)
        virial_avg.append(np.mean(virials))
        virial_err.append(np.std(virials))
        
    ax.errorbar(r_arr,virial_avg,yerr=virial_err,fmt='s',color=colourstr,markersize=5,label=labelstr)

    return r_arr,virial_avg,virial_err



def fit_curve(func,r_arr,sigma_avg,guess,colourstr,labelstr):

    r_val = np.logspace(np.log10(rad[0]*0.8),np.log10(rad[-1]*1.1),100)
    popt,pcov = curve_fit(func,r_arr,virial_avg,p0=guess)
    print(*popt)
    ax.plot(r_val,func(r_val,*popt),color=colourstr,label=labelstr)

    return popt


fig, ax = plt.subplots(ncols=1, sharey=False, subplot_kw=dict(frameon=True), figsize=(5.5,4))


# Initial 
r_arr,virial_avg,virial_err = plot_virial(ax,'virial_term_cloud_20_10_00100.dat','grey','initial')
popt = fit_curve(constant_func,r_arr,virial_avg,(2.46),'pink','initial $2GM/\sigma_{v}^2 R = 4.57\ R^{0}$')
print(*popt)


# After photoion
r_arr,virial_avg,virial_err = plot_virial(ax,'virial_term_cloud_20_10_clsink139_v2_06456.dat','black','after ion')
popt = fit_curve(constant_func,r_arr,virial_avg,(2.35),'red','after ion $2GM/\sigma_{v}^2 R = 4.12\ R^{0}$')
print(*popt)


# Only ionized gas 
r_arr,virial_avg,virial_err = plot_virial(ax,'virial_term_ionized_cloud_20_10_clsink139_v2_06456.dat','navy','ionized gas only')
virial_avg.pop(1)
virial_avg.pop(1)
r_arr.pop(1)
r_arr.pop(1)
popt = fit_curve(powerlaw_func,r_arr,virial_avg,(-1.08,0.084),'royalblue','ionized $2GM/\sigma_{v}^2 R = 0.19\ R^{-0.94}$')
print(*popt)



ax.set_xlabel('size-scale [pc]')
ax.set_ylabel('virial parameter $2GM/\sigma_{v}^2 R$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([1.2e-8,2e3])
ax.legend(loc='lower right')

fig.tight_layout(pad=0.5)
plt.savefig('virial_scale_relation.png',dpi=200)

plt.show()











