# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 19:41:00 2023

@author: lausi
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt



# CMI heating rate 

file = h5py.File("snapshot010.hdf5", "r")
datasetnames = file.keys()
print(file['/Units'].attrs.keys())

udist = file['/Units'].attrs['Unit length in cgs (U_L)']
umass = file['/Units'].attrs['Unit mass in cgs (U_M)']
utime = file['/Units'].attrs['Unit time in cgs (U_t)']

unit_energ = umass*(udist/utime)**2

print('unit_energ',unit_energ)

unit_gamma = unit_energ/utime

coords = np.array(file["/PartType0/Coordinates"])
neutfracH = np.array(file["/PartType0/NeutralFractionH"])

gamma = np.array(file["/PartType0/HeatingRateH"])

file.close()

# cmi ind. cells
nH_cmi_raw = neutfracH
gamma_cmi_raw = gamma*unit_gamma*neutfracH

plt.scatter(nH_cmi_raw,gamma_cmi_raw,s=0.3,color='black')


# bin gamma by nH
ncell = 70
gamma_cells = [0] * ncell
num_cells = [0] * ncell
nH_cells = np.linspace(0,1,ncell)

for ip in range(len(neutfracH)):
    nH = neutfracH[ip]
    icell = np.abs(nH_cells - nH).argmin()
    gamma_cells[icell] += gamma[ip]*unit_gamma * neutfracH[ip]  # to cgs
    num_cells[icell] += 1
    
for icell in range(ncell):
    if num_cells[icell] >= 1:
        gamma_cells[icell] = gamma_cells[icell] / num_cells[icell]
    else:
        gamma_cells[icell] = 0

nH_cells_nozero = []
gamma_cells_nozero = []
for i in range(len(nH_cells)):
    if gamma_cells[i] > 0:
        nH_cells_nozero.append(nH_cells[i])
        gamma_cells_nozero.append(gamma_cells[i])

nH_cmi = nH_cells_nozero
gamma_cmi = gamma_cells_nozero

plt.plot(nH_cmi,gamma_cmi,color='red')


# Phantom heating rate 

data = np.loadtxt('gamma_nH_distri.txt')

nH = data[:,0]
gamma = data[:,1] 

nH_nozero = []
gamma_nozero = []
for i in range(len(nH)):
    if gamma[i] > 0:
        nH_nozero.append(nH[i])
        gamma_nozero.append(gamma[i])

nH_phantom = nH_nozero
gamma_phantom = gamma_nozero
plt.plot(nH_phantom,gamma_phantom,color='blue')


plt.xlabel('nH')
plt.ylabel('gamma [erg/s]')
plt.yscale('log')
plt.show()


fig = plt.figure(figsize=(6, 8),dpi=200)
fig.tight_layout()

ax1 = fig.add_subplot(211)
ax1.set_xlim([0,1])
ax1.set_ylim([1E-25,1E-20])
ax1.set_xlabel('neutral fraction')
ax1.set_ylabel('heating rate [$\mathrm{erg \ s^{-1}}$]')
ax1.set_yscale('log')
ax1.scatter(nH_phantom,gamma_phantom,s=0.5,color='blue',label='Phantom')
ax1.scatter(nH_cmi_raw,gamma_cmi_raw,s=0.3,color='red',label='CMI')
ax1.legend()

ax2 = fig.add_subplot(212)
ax2.set_xlim([0,1])
ax2.set_ylim([1E-25,1E-20])
ax2.set_xlabel('neutral fraction')
ax2.set_ylabel('heating rate [$\mathrm{erg \ s^{-1}}$]')
ax2.set_yscale('log')
ax2.plot(nH_cmi,gamma_cmi,color='red',label='CMI')
ax2.plot(nH_phantom,gamma_phantom,color='blue',label='Phantom')
ax2.legend()


plt.savefig('heating_rate_phantom_cmi')







