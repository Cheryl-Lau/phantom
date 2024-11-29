

import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D


time_cgs = np.loadtxt('nixyzhmf_04999.txt',skiprows=1,max_rows=1)
Myr = 1e6*365*24*60*60 

n,i,x,y,z,h,m,nH = np.loadtxt('nixyzhmf_04999.txt',skiprows=3,unpack=True)
imfilename = 'nixyzhmf_04999.png'


centre = [2.181616076E+00,2.883273141E+00,-2.588239455E+00]  # sink 139
radius = 4.5


fig = plt.figure(figsize=[18,10],dpi=200)
fig.subplots_adjust(wspace=0)


ax1 = fig.add_subplot(122,projection='3d')

ax1.set_ylabel('x [pc]')
ax1.set_xlabel('y [pc]')
ax1.zaxis.set_rotate_label(True)
ax1.set_zlabel('z [pc]', rotation=90)
ax1.view_init(elev=10., azim=330.)

iionized = np.where(nH < 0.8)[0]
x_ionized = x[iionized]
y_ionized = y[iionized]
z_ionized = z[iionized]
nH_ionized = nH[iionized]

pcm = ax1.scatter(x_ionized,y_ionized,z_ionized,c=nH_ionized, cmap='viridis',alpha=0.5, vmin=0, vmax=1)
    
# Get rid of the panes
ax1.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax1.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax1.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

ax1.set_xlim([centre[0]-radius,centre[0]+radius])
ax1.set_xlim([centre[1]-radius,centre[1]+radius])
ax1.set_xlim([centre[2]-radius,centre[2]+radius])


ax2 = fig.add_subplot(121,projection='3d')

ax2.set_ylabel('x [pc]')
ax2.set_xlabel('y [pc]')
ax2.zaxis.set_rotate_label(False)
ax2.set_zlabel('z [pc]', rotation=90)
ax2.view_init(elev=10., azim=240.)

pcm = ax2.scatter(x_ionized,y_ionized,z_ionized,c=nH_ionized, cmap='viridis',alpha=0.5, vmin=0, vmax=1)

ax2.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax2.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax2.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

ax2.set_xlim([centre[0]-radius,centre[0]+radius])
ax2.set_xlim([centre[1]-radius,centre[1]+radius])
ax2.set_xlim([centre[2]-radius,centre[2]+radius])


cbar = fig.colorbar(pcm, ax=[ax1,ax2], location='right', shrink=0.5)
cbar.set_label('neutral fraction', rotation=270, labelpad=20)



    
plt.savefig("ionfrac_particles.png")













