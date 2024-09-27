# -*- coding: utf-8 -*-
"""
Reads rho-gamma-Teq table into pandas 
for extracting information  
"""
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from scipy import interpolate 


rho,gamma,nroot,Teq1,Teq2,Teq3 = np.loadtxt('rho_gamma_Teqs.txt',skiprows=1,unpack=True)

df = pd.DataFrame(columns=('rho','gamma','nroot','Teq1','Teq2','Teq3'))
df['rho'] = rho 
df['gamma'] = gamma 
df['nroot'] = nroot 
df['Teq1'] = Teq1 
df['Teq2'] = Teq2 
df['Teq3'] = Teq3 

# Gamma which corresponds to nH = 0.5 to nH = 0.0 (ionized particles)
gamma_nh00 = 1e-18
gamma_nh05 = 2e-26


#
# Plot regime where there is a 3rd root
#

iselect = np.where((df['nroot']==3) & (df['gamma']>gamma_nh05) & (df['gamma']<gamma_nh00))[0]
rho_select = df['rho'].loc[iselect].to_numpy()
gamma_select = df['gamma'].loc[iselect].to_numpy()
temp3_select = df['Teq3'].loc[iselect].to_numpy()

xmin = np.min(gamma_select)
xmax = np.max(gamma_select)
ymin = np.min(rho_select)
ymax = np.max(rho_select)
zmin = np.min(temp3_select)
zmax = np.max(temp3_select)

#fig,(ax1,ax2) = plt.subplots(nrows=2,sharex=True,subplot_kw=dict(frameon=True),figsize=(8,7),dpi=200)
fig,(ax1,ax2) = plt.subplots(nrows=2,figsize=(7,6),dpi=200)
plt.subplots_adjust(hspace=0.3)

ax1.set_xlabel('Heating rate '+r'$\mathrm{[erg \ s^{-1}]}$')
ax1.set_ylabel('Density '+r'$\mathrm{[g \ cm^{-3}]}$')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim([9e-23,2e-19])
ax1.set_ylim([1e-23,9e-20])
ax1.text(2.8e-21,2e-23,'$\mathrm{2^{nd}}$ stable equilibrium temperature')

x = np.logspace(np.log10(xmin),np.log10(xmax),100)
y = np.logspace(np.log10(ymin),np.log10(ymax),100)
mesh_x,mesh_y = np.meshgrid(x,y)
mesh_z = interpolate.griddata((gamma_select,rho_select), temp3_select, (mesh_x,mesh_y), method='cubic')
mesh_z = np.log10(mesh_z)
pc1 = ax1.pcolormesh(mesh_x,mesh_y,mesh_z,cmap='bwr',vmin=0.8,vmax=9.1)
cbar1 = plt.colorbar(pc1,ax=ax1,pad=0)
cbar1.ax.set_ylabel('log($\mathrm{T_{eq}}$) [K]', rotation=270,labelpad=15.)


#
# Plot regime where they is only 1 root  
#
iselect = np.where((df['nroot']==1) & (df['gamma']>gamma_nh05) & (df['gamma']<gamma_nh00))[0]
rho_select = df['rho'].loc[iselect].to_numpy()
gamma_select = df['gamma'].loc[iselect].to_numpy()
temp1_select = df['Teq1'].loc[iselect].to_numpy()

ax2.set_xlabel('Heating rate '+r'$\mathrm{[erg \ s^{-1}]}$')
ax2.set_ylabel('Density '+r'$\mathrm{[g \ cm^{-3}]}$')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim([9e-23,2e-19])
ax2.set_ylim([1e-23,9e-20])
ax2.text(2.8e-21,2e-23,'$\mathrm{1^{st}}$ stable equilibrium temperature')

x = np.logspace(np.log10(xmin),np.log10(xmax),100)
y = np.logspace(np.log10(ymin),np.log10(ymax),100)
mesh_x,mesh_y = np.meshgrid(x,y)
mesh_z = interpolate.griddata((gamma_select,rho_select), temp1_select, (mesh_x,mesh_y), method='cubic')
mesh_z = np.log10(mesh_z)
pc2 = ax2.pcolormesh(mesh_x,mesh_y,mesh_z,cmap='bwr',vmin=0.8,vmax=9.1) 
cbar2 = plt.colorbar(pc2,ax=ax2,pad=0)
cbar2.ax.set_ylabel('log($\mathrm{T_{eq}}$) [K]', rotation=270,labelpad=15.)


#
# The permitted rho-gamma combo for photoionization heating 
#

def alphaA_func(temp):
        return 6.113e-11 * temp**(-1/2) - 1.857e-13


Tstar = 4e4 
Tstar1 = 9e3
Tstar2 = 6e4

Tstar_range = np.logspace(np.log10(Tstar1),np.log10(Tstar2),1000)

alphaA = alphaA_func(Tstar)
alphaA1 = alphaA_func(Tstar1)
alphaA2 = alphaA_func(Tstar2)

alphaA_range = alphaA_func(Tstar_range)

print(alphaA,alphaA1,alphaA2)

mH = 1.67e-24
kboltz = 1.38e-16

gamma_SN = 0 #1e-21
gamma_bg = 0 #2e-26 

term = (1/mH)*alphaA*(3/2)*kboltz*Tstar
term1 = (1/mH)*alphaA1*(3/2)*kboltz*Tstar1
term2 = (1/mH)*alphaA2*(3/2)*kboltz*Tstar2

term_range = (1/mH)*alphaA_range*(3/2)*kboltz*Tstar_range 

print('term',term,term1,term2)

rho_req = (x -gamma_SN -gamma_bg)/term
rho_req1 = (x -gamma_SN -gamma_bg)/term1
rho_req2 = (x -gamma_SN -gamma_bg)/term2

print('rhos',rho_req[10],rho_req1[10],rho_req2[10])

ax1.plot(x,rho_req,'black',linewidth=0.8,label=r"Permitted $\rho-\Gamma_{photoion}$")
ax2.plot(x,rho_req,'black',linewidth=0.8,label=r"Permitted $\rho-\Gamma_{photoion}$")

term_max = np.max(term_range)
term_min = np.min(term_range)
f1 = lambda x: (x -gamma_SN -gamma_bg)/term_max 
f2 = lambda x: (x -gamma_SN -gamma_bg)/term_min
ax1.fill_between(x, f1(x), f2(x), where=f1(x)<=f2(x), interpolate=True, color='grey', alpha=0.3)
ax2.fill_between(x, f1(x), f2(x), where=f1(x)<=f2(x), interpolate=True, color='grey', alpha=0.3)

ax1.legend(fontsize=9,loc='upper left')
ax2.legend(fontsize=9,loc='upper left')

plt.savefig('temp_roots_HIIregion.png')
plt.show()






















