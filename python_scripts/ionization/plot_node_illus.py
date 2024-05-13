# -*- coding: utf-8 -*-

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors

slice_thickness = 0.05

xyzhmf_parts = np.loadtxt('xyzhmf_allparts.txt',skiprows=3)
x_ap = xyzhmf_parts[:,0]
y_ap = xyzhmf_parts[:,1]
z_ap = xyzhmf_parts[:,2]
h_ap = xyzhmf_parts[:,3]
m_ap = xyzhmf_parts[:,4]
f_ap = xyzhmf_parts[:,5]

x_ap_slice = x_ap[np.where((z_ap>-slice_thickness) & (z_ap<slice_thickness))]
y_ap_slice = y_ap[np.where((z_ap>-slice_thickness) & (z_ap<slice_thickness))]
m_ap_slice = m_ap[np.where((z_ap>-slice_thickness) & (z_ap<slice_thickness))]
f_ap_slice = f_ap[np.where((z_ap>-slice_thickness) & (z_ap<slice_thickness))]


nixyzhmf_beforeiter = np.loadtxt('nixyzhmf_cminode_beforeiterdone.txt',skiprows=1)
x_bi = nixyzhmf_beforeiter[:,2]
y_bi = nixyzhmf_beforeiter[:,3]
z_bi = nixyzhmf_beforeiter[:,4]
h_bi = nixyzhmf_beforeiter[:,5]
m_bi = nixyzhmf_beforeiter[:,6]
f_bi = nixyzhmf_beforeiter[:,7]

x_bi_slice = x_bi[np.where((z_bi>-slice_thickness) & (z_bi<slice_thickness))]
y_bi_slice = y_bi[np.where((z_bi>-slice_thickness) & (z_bi<slice_thickness))]
m_bi_slice = m_bi[np.where((z_bi>-slice_thickness) & (z_bi<slice_thickness))]
f_bi_slice = f_bi[np.where((z_bi>-slice_thickness) & (z_bi<slice_thickness))]


nixyzhmf_afteriter = np.loadtxt('nixyzhmf_cminode_iterdone.txt',skiprows=3)
x_ai = nixyzhmf_afteriter[:,2]
y_ai = nixyzhmf_afteriter[:,3]
z_ai = nixyzhmf_afteriter[:,4]
h_ai = nixyzhmf_afteriter[:,5]
m_ai = nixyzhmf_afteriter[:,6]
f_ai = nixyzhmf_afteriter[:,7]

x_ai_slice = x_ai[np.where((z_ai>-slice_thickness) & (z_ai<slice_thickness))]
y_ai_slice = y_ai[np.where((z_ai>-slice_thickness) & (z_ai<slice_thickness))]
m_ai_slice = m_ai[np.where((z_ai>-slice_thickness) & (z_ai<slice_thickness))]
f_ai_slice = f_ai[np.where((z_ai>-slice_thickness) & (z_ai<slice_thickness))]


fig = plt.figure(figsize=(12, 6.5),dpi=200)
fig.tight_layout()

# masses 

ax1 = fig.add_subplot(231)
ax1.set_xlim([-0.75,0.75])
ax1.set_ylim([-0.75,0.75])
ax1.set_ylabel('y [pc]', fontsize=13)
ax1.scatter(x_ap_slice,y_ap_slice,s=0.3,c=m_ap_slice,cmap='plasma',norm=mcolors.LogNorm())

ax2 = fig.add_subplot(232)
ax2.set_xlim([-0.75,0.75])
ax2.set_ylim([-0.75,0.75])
ax2.scatter(x_bi_slice,y_bi_slice,s=0.3,c=m_bi_slice,cmap='plasma',norm=mcolors.LogNorm())

ax3 = fig.add_subplot(233)
ax3.set_xlim([-0.75,0.75])
ax3.set_ylim([-0.75,0.75])
pcm_m = ax3.scatter(x_ai_slice,y_ai_slice,s=0.3,c=m_ai_slice,cmap='plasma',norm=mcolors.LogNorm())

cbar_m = fig.colorbar(pcm_m, ax=[ax1,ax2,ax3], pad=0.02, location='right')
cbar_m.set_label('Mass [$M_\u2609$]', rotation=270, labelpad=14, fontsize=13)


# neutral frac 

ax4 = fig.add_subplot(234)
ax4.set_xlim([-0.75,0.75])
ax4.set_ylim([-0.75,0.75])
ax4.set_xlabel('x [pc]', fontsize=13)
ax4.set_ylabel('y [pc]', fontsize=13)
ax4.scatter(x_ap_slice,y_ap_slice,s=0.3,c=f_ap_slice,cmap='viridis')

ax5 = fig.add_subplot(235)
ax5.set_xlim([-0.75,0.75])
ax5.set_ylim([-0.75,0.75])
ax5.set_xlabel('x [pc]', fontsize=13)
ax5.scatter(x_bi_slice,y_bi_slice,s=0.3,c=f_bi_slice,cmap='viridis')

ax6 = fig.add_subplot(236)
ax6.set_xlim([-0.75,0.75])
ax6.set_ylim([-0.75,0.75])
ax6.set_xlabel('x [pc]', fontsize=13)
pcm_f = ax6.scatter(x_ai_slice,y_ai_slice,s=0.3,c=f_ai_slice,cmap='viridis')

cbar_f = fig.colorbar(pcm_f, ax=[ax4,ax5,ax6], pad=0.02, location='right')
cbar_f.set_label('Neutral fraction', rotation=270, labelpad=20, fontsize=13)

plt.show()
plt.savefig('node_illus.png')

