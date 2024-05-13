#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot final velocity against s_out for semi-confined case 
"""

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.pyplot import cm
import os 


gamma = 5/3                     # adiabatic index 
G = 6.672041e-8                 # gravitational constant 
kboltz = 1.38066e-16            # Boltzmann constant 
mass_proton = 1.67262158e-24    # proton mass 
gmw = 1.29                      # mean molecular weight 

t_end = 0.8E+14                 # sim end time 
dt = 1E8                        # timestep 
Myr = 1E6*365*24*60*60

r_detect = 25 *3.086e+18        # location of detector 

r_sn = 0.1 *3.086e+18           # SN ejecta radius 
E_sn = 1e51                     # SN energy 
m_sn = 25.0 *1.989e+33          # SN ejecta mass 
vol_sn = 4/3*np.pi*r_sn**3      # volume occupied by SN ejecta 

r_cloud = 24.41 *3.086e+18      # Cloud radius 
rho_cloud = 1e-21               # density of cloud around HII region/cavity
p_cloud = 3.34E-12              # pressure of cloud around HII region/cavity

r_stag = 10. *3.086e+18         # HII region stagnation radius 
#p_hii = p_cloud * 10**(-3*gamma/(1-gamma))          # pressure of HII region from adiabatic relations 
#rho_hii = rho_cloud * (p_hii/p_cloud) * 10**(-3)    # density of HII region from ideal gas law  

rho_inshell = 1.5e-19           # density within swept-up shell 
dr_inshell = 2.0 *3.086e+18     # thickness of shell 

rho_env = 4e-25                 # density of envelope 
p_env = 2.54e-14                # pressure of envelope

omega = 4*np.pi*0.2             # solid angle of channell in Sr



def part_confined_sn(s_out,rho_hii,p_hii):

    # Initial properties of the cavity before being hit by SN shock 
    vol_cav = 4/3*np.pi*r_stag**3
    r_cav = r_stag
    m_inshell = 4*np.pi*r_cav**2*dr_inshell*rho_inshell 
    M_i = 4/3*np.pi*r_cav**3*rho_cloud   # mass of surroundings if it had extended to centre 

    p_cav = (p_hii/(gamma-1)*(vol_cav-vol_sn) + E_sn) * (gamma-1)/vol_cav
    m_cav = m_sn + rho_hii*(vol_cav-vol_sn)
    E_cav = E_sn + p_hii/(gamma-1)*(vol_cav-vol_sn)

    # Initial properties of the shell around the cavity upon hit by SN
    vel_cav = (6*p_cav/rho_inshell)**(1/2)      # velocity by dimension analysis 
    rho_vent = rho_env*(p_cav/p_env)**(1/gamma) # poisson adiabate 
    e_mload = E_cav/m_cav                       # mass loading 

    # Store 
    r_cavi = r_cav 
    m_cavi = m_cav 
    vol_cavi = vol_cav

    t_evol = []
    vout_evol = []

    t = 0.
    while(t < t_end and p_cav > p_env): 
        # Update velocity of escaping gas
        v_out = np.sqrt(2*gamma/(gamma-1)*p_env/rho_env*(((p_cav)/p_env)**((gamma-1)/gamma)-1))  # bernoulli
        if (v_out != v_out):
            print('error sqrt',m_cav,p_cav)
            quit()

        # Update mass contained within cavity 
        dmdt = rho_vent * v_out * s_out 
        m_cav = m_cav - dmdt*dt 
        if (m_cav < 0):
            exit 

        m_shellpswept = m_inshell * (M_i/m_inshell*(r_cav**3/r_cavi**3-1) + 1)           # shell mass + swept-up mass
        dvdt_cav = (1/m_shellpswept) * 4*np.pi*r_cav**2 * (p_cav - rho_cloud*vel_cav**2) # shell EOM  
        dvdt_cav = dvdt_cav - 4/3*np.pi*G*r_cav*rho_vent                                 # take gravity into account
        # Update cavity radius 
        vel_cav = vel_cav + dvdt_cav*dt 
        r_cav = r_cav + vel_cav*dt 

        vol_cav = 4/3*np.pi*r_cav**3 
        
        rho_vent = m_cav/vol_cav 
        p_cav = p_env*(rho_vent/rho_env)**gamma  # poisson adiabate 
        if (p_cav < p_env):
            print('dropped below p_env, stopping')
        ramp_vent = rho_vent * v_out**2
    
        t = t + dt 

        t_evol.append(t)
        vout_evol.append(v_out)

    # later stages 
    icrop = np.where(np.logical_and(np.array(t_evol) > 0.65E+14, np.array(t_evol) < t_end))[0]
    t_evol = np.take(t_evol,icrop)
    vout_evol = np.take(vout_evol,icrop)

    # get mean
    vout = np.mean(vout_evol)
    vouterr = np.std(vout_evol)

    # get gradient 
    coef,cov = np.polyfit(t_evol,vout_evol,1,cov=True)
    grad = coef[0]
    graderr = np.sqrt(cov[0][0])

    # time 
    t_start = t_evol[0]/Myr 
    t_last = t_evol[-1]/Myr 
    print('start and end time [Myr] ',t_start,t_last)

    return vout,vouterr,grad,graderr 

    
def main():

    plt.figure(dpi=300)

    rho_hii_all = np.array([1e-23,1e-22,1e-21,1e-20,1e-19,1e-18])
    #original p_hii for rho_hii = 1e-19:  6.40e-8 

    nomega = 10
    nrho = len(rho_hii_all)
    color = iter(cm.ocean(np.linspace(0.2, 0.8, nrho)))

    # loop over r_detect 
    irun = 1
    for rho_hii in rho_hii_all:

        # calculate corresponding p 
        u = kboltz*1E4/(gmw*mass_proton*(gamma-1))
        p_hii = rho_hii*(gamma-1)*u
        print('> For rho_hii ',rho_hii,'p_hii',p_hii)
        
        filename = 'vout_sout_'+str(int(irun))+'.txt'

        if (os.path.exists(filename) == True):
            omegafrac_test,omega_all,sout_all,vout_all,vouterr_all,voutgrad_all,voutgraderr_all = np.loadtxt(filename,unpack=True)

            c = next(color)
            plt.errorbar(omegafrac_test,vout_all,yerr=vouterr_all,elinewidth=1,color=c,label=str(rho_hii)+' g cm^-3')
    
        else:
            omegafrac_test = np.linspace(0.05,0.6,nomega)
            omega_all = 4*np.pi*omegafrac_test

            # loop over s_out 
            sout_all = []
            vout_all = []
            vouterr_all = []
            voutgrad_all = []
            voutgraderr_all = []

            for omega in omega_all:
                sout = omega*r_stag**2
                print('Running model for omega_frac =',omega/(4*np.pi))
                vout,vouterr,vout_grad,vout_graderr = part_confined_sn(sout,rho_hii,p_hii)
                sout_all.append(sout)
                vout_all.append(vout)
                vouterr_all.append(vouterr)
                voutgrad_all.append(vout_grad)
                voutgraderr_all.append(vout_graderr)
        
        np.savetxt(filename,np.column_stack((omegafrac_test,omega_all,sout_all,vout_all,vouterr_all,voutgrad_all,voutgraderr_all)))
        irun += 1 

    plt.xlabel('Opening [4$\pi$ Sr]')
    plt.ylabel('Outflow velocity [cm s^-1]')
    plt.ylim([0e6,6e6])
#    plt.ylim([3e4,6e6])
#    plt.yscale('log')
    plt.legend()
    plt.savefig('vout_sout_relation.png')
    plt.show()
    
if __name__ == "__main__":
    main()





