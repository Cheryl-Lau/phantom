import numpy as np
import matplotlib.pyplot as plt 
import glob 

partnumsets = ['ten5particles','ten6particles','ten7particles']
totpart = [144666,1426800,9938375]
colours = ['red','blue','black']
labels = [r'$1.47\times10^5$'+' particles', r'$1.43\times10^6$'+' particles', \
            r'$0.99\times10^7$'+' particles']

fig0 = plt.figure()
ax0 = fig0.add_subplot(111)
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
fig7 = plt.figure()
ax7 = fig7.add_subplot(111)
fig8 = plt.figure()
ax8 = fig8.add_subplot(111)



ax0.set_ylabel('total CMI time')
ax0.set_xscale('log')
ax0.set_yscale('log')

ax1.set_ylabel('density func init time')
ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.set_ylabel('Voronoi Density Grid time')
ax2.set_xscale('log')
ax2.set_yscale('log')

ax3.set_ylabel('Grid init time')
ax3.set_xscale('log')
ax3.set_yscale('log')

ax4.set_ylabel('Lloyd time')
ax4.set_xscale('log')
ax4.set_yscale('log')

ax5.set_ylabel('density mapping time')
ax5.set_xscale('log')
ax5.set_yscale('log')

ax6.set_ylabel('sim time')
ax6.set_xscale('log')
ax6.set_yscale('log')

ax7.set_ylabel('photoionization time')
ax7.set_xscale('log')
ax7.set_yscale('log')

ax8.set_ylabel('Reverse mapping time')
ax8.set_xscale('log')
ax8.set_yscale('log')

i = 0

for numpartdir, allpart in zip(partnumsets, totpart):
    path = numpartdir+"/set*"
    for setname in glob.glob(path): 

        # number of pseudo-particles 
        nppartfile = setname+"/nppart"
        nppart = np.loadtxt(nppartfile, unpack=True)

        # cmi runtimes
        cmitimefile = setname+"/ionization-simulation-time-log.txt"
        # entry id	parent id	depth	start time (ticks)	end time (ticks)	start time (s)	end time (s)	label
        start_s,end_s = np.loadtxt(cmitimefile, skiprows=1, unpack=True, usecols=[5,6])

        ax0.scatter(nppart, end_s[-1]-start_s[0], c=colours[i])
        ax1.scatter(nppart, end_s[3]-start_s[3], c=colours[i])
        ax2.scatter(nppart, end_s[4]-start_s[4], c=colours[i])
        ax3.scatter(nppart, end_s[5]-start_s[5], c=colours[i])
        ax4.scatter(nppart, end_s[6]-start_s[6], c=colours[i])
        ax5.scatter(nppart, end_s[7]-start_s[7], c=colours[i])
        ax6.scatter(nppart, end_s[8]-start_s[8], c=colours[i])
        ax7.scatter(nppart, end_s[9]-start_s[9], c=colours[i])
        ax8.scatter(nppart, end_s[10]-start_s[10], c=colours[i])

    i += 1 


plt.show()
