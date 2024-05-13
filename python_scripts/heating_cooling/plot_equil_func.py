# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 


def main():
    mass_proton_cgs = 1.67262158E-24

    temp = np.logspace(2,8,100)
    lambda_cgs = coolcurve(temp)   # erg cm^3 s^-1

    gamma_cgs = 2.5E-19

    rho_cgs = 2.12947615861e-20
    nrho_cgs = rho_cgs/mass_proton_cgs 

    equi_point = gamma_cgs/nrho_cgs
 

    plt.figure(figsize=[7,5],dpi=200)

    plt.plot(temp,lambda_cgs,label='Cooling curve')
    plt.axhline(y = equi_point, color = 'green', linestyle = '-',label='Equilibrium solution')
  
    plt.xlabel('T [K]')
    plt.ylabel("$\Lambda$ [$ \mathrm{ erg \ cm^{3} \ s^{-1} } $] ")
    
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1E-24,1E-20])
    plt.legend()
    plt.show()
    
    plt.savefig('roots_on_coolcurve.png')
    
    




def coolfunc_mod_derijcke(temp): 
    '''
    Modified JML06 cooling function to mimic DR13
    [in cgs units]
    '''
    GammaKI_cgs = 2E-26  

    # First term 
    if temp > 10**4.15: 
        lambdagamma1 = 4.69414E-4 * 1E7 * np.exp(-1.184E5 * 1.15983E6 / ((temp*10**(-0.08))**2.68935 + 1000)) 
    else: 
        lambdagamma1 = 1E7 * np.exp(-1.184E5 / ((temp*10**(-0.16)) + 1000))* temp**0.18 
    
    # Second term
    lambdagamma2 = 0.115 * 0.014 * (temp*10**(-0.75))**0.60 * np.exp(-72/(temp*10**(-0.12))) #*10**0.2
    
    lambdagamma2 = lambdagamma2 * temp**(0.15)
    
    # shift up 
    lambdagamma = (lambdagamma1 + lambdagamma2) * 10**(0.95)
    
    # correct low temp part 
    if temp < 10**3.7:
        lambdagamma = lambdagamma * 10**1.665 * temp**(-0.45)
    
    # High temp part
    if temp > 10**5.3: 
        if temp < 10**6.5: 
            lambdagamma = lambdagamma * 10**5.3 * temp**(-1.)
        else: 
            lambdagamma = lambdagamma * 10**0.425 * temp**(-0.25) 
            
    # Highest temp part
    if temp > 10**7.5:
        lambdagamma = lambdagamma * 10**(-5.25) * temp**(0.7)
     
    # shift down 
    lambdagamma = lambdagamma / 10**(0.35)
        
     
    lambdacool = lambdagamma * GammaKI_cgs   
    
    return lambdacool 


def coolcurve(temp):
    lambda_whole = []
    for i in range(len(temp)): 
        lambda_pt = coolfunc_mod_derijcke(temp[i])
        lambda_whole.append(lambda_pt)
    return lambda_whole

main()