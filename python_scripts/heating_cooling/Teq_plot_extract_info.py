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
gamma_nh00 = 2.8e-21 
gamma_nh05 = 7.9e-22 

iselect = np.where((df['nroot']==3) & (df['gamma']>gamma_nh05) & (df['gamma']<gamma_nh00))[0]


rho_select = df['rho'].loc[iselect].to_numpy()
gamma_select = df['gamma'].loc[iselect].to_numpy()


#
# Plot 2nd root - going above which will drift to 3rd root 
#
temp2_select = df['Teq2'].loc[iselect].to_numpy()

xmin = np.min(gamma_select)
xmax = np.max(gamma_select)
ymin = np.min(rho_select)
ymax = np.max(rho_select)
zmin = np.min(temp2_select)
zmax = np.max(temp2_select)

fig,(ax1,ax2) = plt.subplots(nrows=2,sharex=True,subplot_kw=dict(frameon=True),figsize=(8,10),dpi=200)
plt.subplots_adjust(hspace=.0)

ax1.set_xlabel('Heating rate [erg/s]')
ax1.set_ylabel('Density [g cm^-3]')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim([7e-22,3.1e-21])
ax1.set_ylim([9e-24,4e-22])
ax1.text(8e-22,3e-22,'Threshold temperature')

x = np.logspace(np.log10(xmin),np.log10(xmax),100)
y = np.logspace(np.log10(ymin),np.log10(ymax),100)
mesh_x,mesh_y = np.meshgrid(x,y)
mesh_z = interpolate.griddata((gamma_select,rho_select), temp2_select, (mesh_x,mesh_y), method='cubic')
mesh_z = np.log10(mesh_z)
pc1 = ax1.pcolormesh(mesh_x,mesh_y,mesh_z,cmap='viridis',vmin=5.8,vmax=9.1)
#cbar1 = plt.colorbar(pc1,ax=ax1)
#cbar1.ax.set_ylabel('log(T) [K]', rotation=270,labelpad=15.)



#
# Plot 3rd root - the final temp that it will reach 
#
temp3_select = df['Teq3'].loc[iselect].to_numpy()

ax2.set_xlabel('Heating rate [erg/s]')
ax2.set_ylabel('Density [g cm^-3]')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim([7e-22,3.1e-21])
ax2.set_ylim([9e-24,4e-22])
ax2.text(8e-22,3e-22,'2nd stable equilibrium temperature')

x = np.logspace(np.log10(xmin),np.log10(xmax),100)
y = np.logspace(np.log10(ymin),np.log10(ymax),100)
mesh_x,mesh_y = np.meshgrid(x,y)
mesh_z = interpolate.griddata((gamma_select,rho_select), temp3_select, (mesh_x,mesh_y), method='cubic')
mesh_z = np.log10(mesh_z)
pc2 = ax2.pcolormesh(mesh_x,mesh_y,mesh_z,cmap='viridis',vmin=5.8,vmax=9.1) 
cbar2 = plt.colorbar(pc2,ax=(ax1,ax2),pad=0)
cbar2.ax.set_ylabel('log(T) [K]', rotation=270,labelpad=15.)


plt.savefig('temp_roots_HIIregion.png')
plt.show()






















