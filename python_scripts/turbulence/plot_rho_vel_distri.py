
import numpy as np 
import matplotlib.pyplot as plt 

rho_cgs,num_rho = np.loadtxt('density_binned_cloud_20_10_clsink139_v2_06456.dat',unpack=True)

vel_cgs,num_vel = np.loadtxt('velocity_binned_cloud_20_10_clsink139_v2_06456.dat',unpack=True)




fig, (ax1,ax2) = plt.subplots(ncols=2, sharey=False, subplot_kw=dict(frameon=True), figsize=(10,5))

ax1.plot(rho_cgs,num_rho)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('density [$\mathrm{g\ cm^{-3}}$]')
ax1.set_ylabel('number of particles')

ax2.plot(vel_cgs,num_vel)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('velocity [$\mathrm{cm\ s^{-1}}$]')
ax2.set_ylabel('number of particles')

fig.tight_layout(pad=0.5)
plt.savefig('vel_rho_distri.png',dpi=200)


plt.show()
