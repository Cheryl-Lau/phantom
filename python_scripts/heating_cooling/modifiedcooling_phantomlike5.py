# -*- coding: utf-8 -*-


import numpy as np 
import pandas as pd 
import random 


# Constants from phantom 
GammaKI_cgs = 2E-26  
mass_proton_cgs = 1.67262158E-24 
rhominJML_cgs = 7E-29  # min rho of which there're roots 
rhomaxJML_cgs = 1E-19
maxt = 1000
kboltz = 1.38066E-16
gmw = 2.381 
gamma = 5./3.

# Code units (setup_sphereinbox)
unit_density = 6.768E-23
unit_velocity = 6.558E3
utime = 4.706E14
udist = 3.086E18
umass = 1.989E33

# Range in temp from Joung & Mac Low 
TminJML = 1E1
TmaxJML = 1E8
    
# Init tables 
nteq_table = np.zeros((maxt,5), dtype=object)
rhoueq_table = np.zeros((maxt,5), dtype=object)    

# Cooling params in code units (consistent with phantom) 
LambdaKI_coef = GammaKI_cgs*umass*utime**3/(mass_proton_cgs**2*udist**5) 
GammaKI = GammaKI_cgs*utime**3/(mass_proton_cgs*udist**2)   


# Test with real particle properties  
rho_parts = []
temp_parts = [] 
tau_parts = [] 
eni_parts = [] 
eninew_parts = [] 



# Subroutine 
def init_hv4tableJML(): 
    '''
    Corresponds to init_hv4table in phantom 
    Creates a table of nrho vs equitemps 
    Then convert to rho vs equius 
    '''

    nrhomin_cgs = rhominJML_cgs/mass_proton_cgs 
    nrhomax_cgs = rhomaxJML_cgs/mass_proton_cgs 
    dnrho_cgs = (np.log10(nrhomax_cgs)-np.log10(nrhomin_cgs))/maxt 
    
    # Local max and min of cooling curve
    T01 = 198609.4917357372
    T02 = 31622776.6016837933
    
    
    # Solve for Teq1, Teq2 and Teq3 for each value of n 
    for i in range(maxt): 
        nrho_cgs = nrhomin_cgs * 10**(i*dnrho_cgs)  # fortran: (i-1)*dnrho_cgs
        Teq1,irterr1 = root_bisection(nrho_cgs,TminJML,T01)
        Teq2,irterr2 = root_bisection(nrho_cgs,T01,T02)
        Teq3,irterr3 = root_bisection(nrho_cgs,T02,TmaxJML) 
        
        if irterr3 == 1:
            if irterr2 == 1:
                if irterr1 == 1:
                    numroots = 0
                else:
                    numroots = 1
            else:
                numroots = 2
        else:
            numroots = 3
        
        nteq_table[i][0] = nrho_cgs 
        nteq_table[i][1] = numroots 
        nteq_table[i][2] = Teq1 
        nteq_table[i][3] = Teq2 
        nteq_table[i][4] = Teq3
        
        # Convert to rho and u in code units    
        rho = nrho_cgs * mass_proton_cgs / unit_density 
        ueq1 = kboltz * Teq1 / (gmw*mass_proton_cgs*(gamma-1.)) / unit_velocity**2  
        ueq2 = kboltz * Teq2 / (gmw*mass_proton_cgs*(gamma-1.)) / unit_velocity**2
        ueq3 = kboltz * Teq3 / (gmw*mass_proton_cgs*(gamma-1.)) / unit_velocity**2
        rhoueq_table[i][0] = rho 
        rhoueq_table[i][1] = numroots 
        rhoueq_table[i][2] = ueq1 
        rhoueq_table[i][3] = ueq2  
        rhoueq_table[i][4] = ueq3  
    
    print(pd.DataFrame(nteq_table,columns=['n','numroots','Teq1','Teq2','Teq3']))
    print(pd.DataFrame(rhoueq_table,columns=['rho','numroots','ueq1','ueq2','ueq3']))

    with open('rhoueq_table.txt', 'w') as f:
        pd.set_option('display.max_rows',None)
        f.write(str(pd.DataFrame(rhoueq_table,columns=['rho','numroots','ueq1','ueq2','ueq3'])))
        pd.reset_option('display.max_rows')
        
    return 
    


# Real function 
def coolfunc(temp): 
    '''
    Modified VS07 cooling function to mimic JML06  
    [in cgs units]
    '''

    # First term 
    if temp > 10**4.15: 
        lambdagamma1 = 4.69414E-4 * 1E7 * np.exp(-1.184E5 * 1.15983E6 / ((temp*10**(-0.08))**2.68935 + 1000)) 
    else: 
        lambdagamma1 = 1E7 * np.exp(-1.184E5 / ((temp*10**(-0.16)) + 1000))* temp**0.18 
    
    # Second term
    lambdagamma2 = 0.215 * 0.014 * (temp*10**(-0.2))**0.66 * np.exp(-92/(temp*10**(-0.12))) *10**0.2
    
    # shift up 
    lambdagamma = (lambdagamma1 + lambdagamma2) * 10**(0.75)
    
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
     
    lambdacool = lambdagamma * GammaKI_cgs   
    
    return lambdacool  


    


# Real function 
def equifunc(nrho,temp): 
    '''
    f(T) = nrho*lambda(T) - gamma  [in cgs units] 
    '''
    lambdacool = coolfunc(temp)  
    f_T = nrho*lambdacool - GammaKI_cgs
    
    return f_T



# Subroutine 
def root_bisection(nrho,Tmin,Tmax):
    '''
    Brackets Teq with Tmin and Tmax 
    Returns irterr = 1 and Teq = 0 if no roots found
    '''
    n_iter = 0
    tol = 2
    converged = False 
    irterr = 0   # 1 indicates no roots found 
    
    # Check signs of input Tmin and Tmax 
    Tminsign = np.sign(equifunc(nrho,Tmin))  # fortran: sign(1,equifunc(nrho,Tmin)) 
    Tmaxsign = np.sign(equifunc(nrho,Tmax))
    
    if Tminsign == Tmaxsign: 
        Teq = 0
        irterr = 1 

    while converged==False and irterr == 0: 
        
        if (Tmax-Tmin) < tol:
            converged = True 
            Teq = (Tmax+Tmin)/2. 
            
        Tmid = (Tmax+Tmin)/2. 
        f_T = equifunc(nrho,Tmid) 
        
        if Tminsign == -1:    # ascending
            if f_T >= 0:
                Tmax = Tmid 
            elif f_T < 0:
                Tmin = Tmid 
        elif Tmaxsign == -1:  # descending 
            if f_T < 0:
                Tmax = Tmid 
            elif f_T >= 0:
                Tmin = Tmid 
        n_iter += 1 
        
        if n_iter > 100: 
            Teq = 0.
            irterr = 1 
            raise Exception('something wrong in root-search')
            
    return Teq, irterr 



# Subroutine 
def joungmaclow_implicit(rhoi,eni): 
    '''
    Corresponds to koyamaintsuka implicit cooling subroutine in phantom 
    With the input rho and internal energy of particle, 
    interpolate from table to compute its ueq 
    '''
    
    rhotable = rhoueq_table[:,0] 
    numroottable = rhoueq_table[:,1]
    
    
    if rhoi < rhotable[0]:  # fortran: [1] 
        # diectly use fist entry 
        numroots = numroottable[0] 
        ueqs = [rhoueq_table[0,2], rhoueq_table[0,3], rhoueq_table[0,4]] 
        
    else: 
    
        # minloc of input rho
        i = np.argmin(abs(np.array(rhotable)-rhoi))  
        
        # [j] lower bound; [j+1] upper bound 
        
        # Avoid out of bounds 
        if i == 0:          # fortran: i == 1
            j = 0           # fortran: j = 1 
        elif i == maxt-1:   # fortran: i == maxt 
            # extrapolate using the last two entries 
            j = i-1         # fortran: j = i-1   
            
        elif rhotable[i] >= rhoi and rhotable[i-1] < rhoi: 
            j = i-1 
        elif rhotable[i] < rhoi and rhotable[i+1] >= rhoi: 
            j = i 
        
        
        ### Interpolate to give the final ueqs ### 
        # Note: Interpolate only if both ueq[j] and ueq[j+1] exist  
        #       otherwise, use the nearest root. 
        
        numroots = numroottable[i]  
        if numroots == 0:
            raise ValueError('no roots found for one particle')
        
        # Init to store the interpolated ueqs 
        ueqs = np.zeros((numroots))
        
        for r in range(1,numroots+1): 
                 
            utable = rhoueq_table[:,1+r]  # fortran: [2+r] 
            
            if utable[j] > 0. and utable[j+1] > 0.:            
                
                ueqs[r-1] = utable[j] + (rhoi-rhotable[j])*(utable[j+1] - utable[j])/(rhotable[j+1]-rhotable[j])  
                # fortran: ueqs[r]
                  
            else: 
    
                if utable[j] == 0. and utable[j+1] > 0.: 
                    ueqs[r-1] = utable[j+1]  
                elif utable[j] > 0. and utable[j+1] == 0.: 
                    ueqs[r-1] = utable[j] 
                    # fortran: ueqs[r]
                else: 
                    raise Exception('unexpected cases exist')
                
         
    
    ### Setting the right ueq with the input eni ###
    if numroots == 1: 
        ueq_final = ueqs[0] 
        
    else: 
        if eni <= ueqs[1]: # fortran: ueqs[2] 
            ueq_final = ueqs[0]  # fortran: ueqs[1]
        else:
            if numroots == 2:
                ueq_final = kboltz*TmaxJML / (gmw*mass_proton_cgs*(gamma-1)) / unit_velocity**2
            elif numroots == 3: 
                ueq_final = ueqs[2]  # fortran: ueqs[3] 
            
    
    # Isothermal temp of particle
    tempi = gmw*mass_proton_cgs/kboltz*(gamma-1.)*eni * unit_velocity**2 
    
    
    # Test timescales 
    rho_parts.append(rhoi)
    temp_parts.append(tempi)
    

    ### Compute new internal energy ###  
    LambdaKI = LambdaKI_coef * coolfunc(tempi)/GammaKI_cgs
    dudti = rhoi*LambdaKI - GammaKI 
    deni = eni - ueq_final 
    
    dt = 1E-5 
    
    if abs(deni) > 0.: 
        tau = abs(deni/dudti)  
        eni_new = ueq_final + deni*np.exp(-dt/tau) 
        dudti = -(eni-eni_new)/dt  
        
        eni_parts.append(eni) 
        eninew_parts.append(eni_new)
        tau_parts.append(tau) 
        
    else: # no change in u 
        dudti = -dudti 

    return dudti 




def main(): 
    
    
    test = 2
        
    
    # Phantom init stage
    init_hv4tableJML() 
    
    
    # Phantom runtime stage 
    if test == 1: 
        
        rho = 1000.  # 0.00019 # 0.00001  #6.1 
        en = 100.
        joungmaclow_implicit(rho,en) 
        
        
    elif test == 2: 
        
        # Particle properties [code units] 
        en = 1E5
        rhos = np.logspace(-1,3,100) 
        for rho in rhos: 
            joungmaclow_implicit(rho,en) 
            
            
        rhocgs_parts = [i*unit_density for i in rho_parts]  
        nrhocgs_parts = [i/mass_proton_cgs for i in rhocgs_parts]
        enicgs_parts = [i*unit_velocity**2 for i in eni_parts] 
        eninewcgs_parts = [i*unit_velocity**2 for i in eninew_parts] 
        
        pd.options.display.max_columns = 10
        part_table = zip(temp_parts,rhocgs_parts,enicgs_parts,eninewcgs_parts,tau_parts)
        print(pd.DataFrame(part_table,columns=['T [K]','rho [g/cm^3]','old u [erg/g]','new u [erg/g]','timescale']))
        
        
    elif test == 3: 
        
        numpart = 10 
        for ip in range(numpart): 
        
            # Particle properties [code units] from sphngsne_test13  
            rho = random.uniform(0.,10**3) 
            en = random.uniform(10.,10**3.5)  
            joungmaclow_implicit(rho,en) 
        
    return 



if __name__ == "__main__":
    main()
    

    
    
    
    
    
    
    

