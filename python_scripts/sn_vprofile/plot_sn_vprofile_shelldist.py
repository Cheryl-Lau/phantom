
import numpy as np 
import matplotlib.pyplot as plt 

Myr = 1E6*365*24*60*60 
pc = 3.086e+18


def main(): 

    rad, vel = np.loadtxt('sn_velocity_profile.dat',skiprows=1,unpack=True)  # cgs units 

    tmax = 0.0007*Myr
    taxis = np.linspace(0,tmax,100)
    
    fig,ax = plt.subplots(nrows=1, sharey=True, subplot_kw=dict(frameon=True), figsize=(5.5,4), dpi=150)

    n = 0
    icount = 0
    for i in range(len(rad)):
        for j in range(len(rad)):   
            idone = n/(len(rad))**2*100
            if (idone > icount):
                icount += 10
                print(int(idone),'% done')
            r1 = rad[i]
            r2 = rad[j]
            v1 = vel[i]
            v2 = vel[j]
            if ((v1 > v2) and (r1 < r2)): 
                dr = abs(rad[i] - rad[j])
                r_shell = dr*v1/(v1-v2)  # dist to which v will hit
                t_shell = r_shell/v1     # time of hitting 
                ax.scatter(r_shell/pc,t_shell/Myr,s=2,color='red')
                ax.plot(trajectory(r_shell,t_shell,(v1+v2)/2,taxis)/pc,taxis/Myr,linewidth=0.1,color='coral',alpha=0.1)
            elif ((v1 > v2) and (r1 > r2)): 
                ax.plot(trajectory(r1,0,v1,taxis)/pc,taxis/Myr,linewidth=0.1,color='grey',alpha=0.5)
                ax.plot(trajectory(r2,0,v2,taxis)/pc,taxis/Myr,linewidth=0.1,color='grey',alpha=0.5)
            n += 1 

    ax.set_ylabel('time [Myr]')
    ax.set_xlabel('distance [pc]')
    ax.set_xlim([0,1.5])
    fig.tight_layout(pad=0.5)

    plt.savefig('sn_vprofile_shells.png')

    plt.show()

    return 


def trajectory(r0,t0,v,t):
    return v*(t-t0) + r0 


main()
    
