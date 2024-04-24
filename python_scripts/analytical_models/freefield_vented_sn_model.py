#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare the time evolution of outflow velocity and pressure for 
free-field supernova and partially-confined supernova
"""

import numpy as np 
import matplotlib.pyplot as plt 

plot_velocity = True 
plot_thermp = True 
plot_ramp = True 

# Sedov is a self-similar solution, here we'll put everything in SI units 

gamma = 5/3                     # adiabatic index 
G = 6.67e-11                    # gravitational constant 

t_end = 1.5E13                  # sim end time 
dt = 1E8                        # timestep 

r_detect = 30. *3.086e+16       # location of detector 

r_sn = 0.1 *3.086e+16           # SN ejecta radius 
E_sn = 0.7*1e51 *1e-7           # SN energy 
m_sn = 1e8 *1.989e+30           # SN ejecta mass 
vol_sn = 4/3*np.pi*r_sn**3      # volume occupied by SN ejecta 

r_cloud = 24.41 *3.086e+16      # Cloud radius 
rho_cloud = 1e-21 *1e+3         # density of cloud around HII region/cavity
p_cloud = 6.36e-13 *1e-1        # pressure of cloud around HII region/cavity

r_stag = 10. *3.086e+16         # HII region stagnation radius 
p_hii = p_cloud * 10**(-3*gamma/(1-gamma))          # pressure of HII region from adiabatic relations 
rho_hii = rho_cloud * (p_cloud/p_hii)**(1/gamma)    # density of HII region from adiabatic relations 

rho_inshell = 1.5e-19           # density within swept-up shell 
dr_inshell = 2.0 *3.086e+16     # thickness of shell 

rho_env = 4e-25 *1e+3           # density of envelope 
p_env = 2.54e-14 *1e-1          # pressure of envelope 

omega = 4*np.pi*0.1             # solid angle of channell in Sr

# NOTE: p - thermal pressure; ramp - ram pressure 
#       if assuming SN energy is completely thermalized by the time it reaches the shell,
#       then will only need to consider thermal pressure p


'''
 - Free-field explosion -
 A SN injected into the envelope (low-density ISM) of rho_env, p_env.
 We place a detector at distance r_detect.
 Suppose SN shock shell already reached r_detect. 
'''
def free_field_sn(ax1,ax2,ax3):

    place_in_HIIregion = False
    place_in_cloud = True 

    s_out = omega*r_detect**2   # to compare with confined case 
    
    time = []
    v_detect_evol = []          # velocity at detector
    p_detect_evol = []          # thermal pressure 
    ramp_detect_evol = []       # ram pressure 
    
    t = 0.
    etot = 0.                   # accumulated energy passing through s_out 
    r_shell = r_detect          # initial radius of shell 

    while (t < t_end): 

        global rho_env, p_env

        if (place_in_HIIregion == True):
            if (r_shell < r_stag):
                rho_env = rho_hii
                p_env = p_hii
                print('in HII') 
            elif (r_shell < r_stag+dr_inshell): 
                rho_env = rho_inshell 
                p_env = p_inshell
                print('in shell')
            elif (r_shell < r_cloud):
                rho_env = rho_cloud
                p_env = p_cloud 
        elif (place_in_cloud == True):
            if (r_shell < r_cloud):
                rho_env = rho_cloud 
                p_env = p_cloud 

        # Update shell radius as it expands 
        v_shell = 2/5 * (1.17)**(5/2) * (E_sn/rho_env)**(1/2) * r_shell**(-3/2)  # sedov sol 
        r_shell = r_shell + v_shell*dt 
        
        # Properties immediately behind shock according to Sedov equations
        v_shell = 2/5 * (1.17)**(5/2) * (E_sn/rho_env)**(1/2) * r_shell**(-3/2)
        p_shell = 4/25 * (1.17)**5 * E_sn * r_shell**(-3)
        rho_shell = (1.17/r_shell)**5 * E_sn * t**2
        cs = np.sqrt(gamma*p_env/rho_env)
        mach = v_shell/cs
        
        # Properties within shocked region 
        r_fac = r_detect/r_shell
        a = (7*gamma-1)/(gamma**2-1)
        b = 3/(gamma-1)
        c = (2*gamma+10)/(gamma-7)
        d = (2*gamma**2+7*gamma-3)/(gamma-7)
        v_detect = v_shell * (r_fac/gamma + ((gamma-1)/(gamma**2+gamma)) * r_fac**a)
        rho_detect = rho_env * (gamma+1)/(gamma-1)*r_fac**b/gamma**c*(gamma+1-r_fac**(a-1))**c
        p_detect = p_env * mach**2 * (2*gamma**(1-d))/(gamma+1) * (gamma+1-r_fac**(a-1))**d
        ramp_detect = rho_detect * v_detect**2
        
        # Energy passing through area s_out at r_detect within dt 
        m_box = rho_detect * s_out * v_detect * dt 
        e_box = 0.5*m_box*v_detect**2 + p_detect/(rho_detect*(gamma-1))*m_box
        etot = etot + e_box 
        
        t = t + dt 
        time.append(t)
        v_detect_evol.append(v_detect)
        p_detect_evol.append(p_detect)
        ramp_detect_evol.append(ramp_detect)
        
        
    if (plot_velocity):
        ax1.scatter(time,v_detect_evol,s=0.1,color='blue',label='free-field')
        ax1.set_yscale('log')
        ax1.set_xlabel('time [s]')
        ax1.set_ylabel('velocity [m/s]')
        
    if (plot_thermp):
        ax2.scatter(time,p_detect_evol,s=0.1,color='blue',label='free-field')
        ax2.set_yscale('log')
        ax2.set_xlabel('time [s]')
        ax2.set_ylabel('thermal pressure [kg/m/s^2]')

    if (plot_ramp):
        ax3.scatter(time,ramp_detect_evol,s=0.1,color='blue',label='free-field')
        ax3.set_yscale('log')
        ax3.set_xlabel('time [s]')
        ax3.set_ylabel('ram pressure [kg/m/s^2]')

    print('- Free-field SN -')
    print('Supernova energy: ',E_sn)
    print('Supernova energy in direction of detector: ',E_sn*omega/(4*np.pi))
    print('Accumulated energy in direction of detector: ',etot,' J')



'''
 - Partially-confined explosion -
 A SN injected in a cavity within a uniform spherical molecular cloud (rho_cloud) carved by 
 photoionization, surrounded by a dense shell of gas. 
 The cloud is enveloped by low-density hot ISM of rho_env, p_env.
 The cavity further expands upon being hit by SN shock and sweeps up more gas into the shell. 
 A channel of cross-section area s_out bridges the cavity and the envelope, allowing gas to 
 escape, passing through the cloud. 
 We assume the gas properties throughout this vent remain the same - neglecting KH instability.
 We place a detector at the exit of the vent r_detect = r_cloud
 Suppose SN shock is already completely thermalized within the cavity.  
'''
def part_confined_sn(with_HII,expand_cav,ax1,ax2,ax3):
    
    s_out = omega*r_detect**2

    # Initial properties of the cavity before being hit by SN shock 
    vol_cav = 4/3*np.pi*r_stag**3
    r_cav = r_stag
    m_inshell = 4*np.pi*r_cav**2*dr_inshell*rho_inshell 
    M_i = 4/3*np.pi*r_cav**3*rho_cloud   # mass of surroundings if it had extended to centre 
    
    if (with_HII):
        p_cav = (p_hii/(gamma-1)*(vol_cav-vol_sn) + E_sn) * (gamma-1)/vol_cav
        m_cav = m_sn + rho_hii*(vol_cav-vol_sn)
        E_cav = E_sn + p_hii/(gamma-1)*(vol_cav-vol_sn)
    else: 
        p_cav = E_sn * (gamma-1)/vol_cav
        m_cav = m_sn + rho_env*(vol_cav-vol_sn)
        E_cav = E_sn
    
    # Initial properties of the shell around the cavity upon hit by SN
    vel_cav = (6*p_cav/rho_inshell)**(1/2)      # velocity by dimension analysis 
    rho_vent = rho_env*(p_cav/p_env)**(1/gamma) # poisson adiabate 
    e_mload = E_cav/m_cav                       # mass loading 

    # Store 
    r_cavi = r_cav 
    m_cavi = m_cav 
    vol_cavi = vol_cav

    # Init ram pressure at vent to aid the escape 
    v_shell = 2/5 * (1.17)**(5/2) * (E_sn/rho_env)**(1/2) * r_cav**(-3/2)
    ramp_vent = rho_vent * v_shell**2

    time = []
    v_vent_evol = []        # velocity at vent / detector
    p_vent_evol = []        # thermal pressure 
    ramp_vent_evol = []     # ram pressure 
    vel_cav_evol = []       # velocity of cavity shell kicked by SN shock
    r_cav_evol = []         # radius of cavity

    t = 0.
    etot = 0.               # accumulated energy passing through s_out 
    mtot = 0.               # accumulated mass

    while(t < t_end): 
        # Update velocity of escaping gas
        v_out = np.sqrt(2*gamma/(gamma-1)*p_env/rho_env*(((p_cav)/p_env)**((gamma-1)/gamma)-1))  # bernoulli

        # Update mass contained within cavity 
        dmdt = rho_vent * v_out * s_out 
        m_cav = m_cav - dmdt*dt 
        
        if (expand_cav):
            m_shellpswept = m_inshell * (M_i/m_inshell*(r_cav**3/r_cavi**3-1) + 1)           # shell mass + swept-up mass
            dvdt_cav = (1/m_shellpswept) * 4*np.pi*r_cav**2 * (p_cav - rho_cloud*vel_cav**2) # shell EOM
            dvdt_cav = dvdt_cav - 4/3*np.pi*G*r_cav*rho_vent                                 # take gravity into account
            # Update cavity radius 
            vel_cav = vel_cav + dvdt_cav*dt 
            r_cav = r_cav + vel_cav*dt 
        
        vol_cav = 4/3*np.pi*r_cav**3 
        
        rho_vent = m_cav/vol_cav 
        p_cav = p_env*(rho_vent/rho_env)**gamma  # poisson adiabate 
        ramp_vent = rho_vent * v_out**2
    
        # Energy passing through the vent of area s_out within dt
        m_box = rho_vent * s_out * v_out * dt 
        e_box = m_box * e_mload 
        etot = etot + e_box 
        mtot = mtot + m_box 

        t = t + dt 
        time.append(t)
        v_vent_evol.append(v_out)
        p_vent_evol.append(p_cav)
        ramp_vent_evol.append(ramp_vent)
        if (expand_cav):
            vel_cav_evol.append(vel_cav)
            r_cav_evol.append(r_cav)
        
    
    if (plot_velocity):
        ax1.scatter(time,v_vent_evol,s=0.1,color='red',label='partially-confined')
        ax1.set_yscale('log')
        ax1.set_xlabel('time [s]')
        ax1.set_ylabel('velocity [m/s]')   

    if (plot_thermp):
        ax2.scatter(time,p_vent_evol,s=0.1,color='red',label='partially-confined')
        ax2.set_yscale('log')
        ax2.set_xlabel('time [s]')
        ax2.set_ylabel('thermal pressure [kg/m/s^2]')  
        
    if (plot_ramp):
        ax3.scatter(time,ramp_vent_evol,s=0.1,color='red',label='partially-confined')
        ax3.set_yscale('log')
        ax3.set_xlabel('time [s]')
        ax3.set_ylabel('ram pressure [kg/m/s^2]')
    
    if (expand_cav):
        fig4 = plt.figure(figsize=[7,5])
        ax4 = fig4.add_subplot(111)
        ax4.scatter(time,r_cav_evol,s=0.1)
        ax4.set_xlabel('time [s]')
        ax4.set_ylabel('cavity shell radius [m]')

        fig5 = plt.figure(figsize=[7,5])
        ax5 = fig5.add_subplot(111)
        ax5.scatter(time,vel_cav_evol,s=0.1)
        ax5.set_yscale('log')
        ax5.set_xlabel('time [s]')
        ax5.set_ylabel('cavity shell velocity [m/s]')

    print('')
    print('- Partially-confined SN -')
    print('Supernova energy: ',E_sn)
    if (with_HII):
        print('HII region energy: ',p_hii/(gamma-1)*(vol_cavi-vol_sn))
    print('Accumulated escaped energy: ',etot,' J')
    print('Initial mass: ',m_cavi,' kg')
    print('Accumulated mass: ',mtot,' kg')    
    print('Initial cavity radius',r_cavi,'m')
    print('Final cavity radius',r_cav,'m')
    
    
def main():
     
    with_HII_region = True 
    expanding_cavity = True

    fig1 = plt.figure(figsize=[7,5])    # vx
    ax1 = fig1.add_subplot(111)
#    ax1.set_ylim([5E3,6E5])

    fig2 = plt.figure(figsize=[7,5])    # thermp
    ax2 = fig2.add_subplot(111)
#    ax2.set_ylim([2E-13,6E-6])

    fig3 = plt.figure(figsize=[7,5])    # ramp
    ax3 = fig3.add_subplot(111)
#    ax3.set_ylim([3E-15,5E-5])
    
    free_field_sn(ax1,ax2,ax3)
    part_confined_sn(with_HII_region,expanding_cavity,ax1,ax2,ax3)
    
    ax1.legend()
    ax2.legend()
    ax3.legend()

    fig1.savefig('mod_vx.png', dpi=200)
    fig2.savefig('mod_thermp.png', dpi=200)
    fig3.savefig('mod_ramp.png', dpi=200)

    plt.show()
    
    
if __name__ == "__main__":
    main()
        





