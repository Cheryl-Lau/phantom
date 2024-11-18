# -*- coding: utf-8 -*-


import numpy as np 
import random as rd 
import pandas as pd 


'''
Construct table of a, m, Nsn, Ek which satisfies the 
relation sum_i_Nsn 1/2mv^2 = Ek 
'''

unit_energ = 8.554E+40 
umass = 1.989E+33

test_one = False

# Params to scan 
a_scan = np.linspace(2E9,3E9,200)
Nsn_scan = np.linspace(10,3000,200)
mi_scan = np.linspace(1E-3,1E-1,200)


if test_one == True: 
    a_scan = [2E9,3E9]
    Nsn_scan = [10,3000]
    mi_scan = [1E-3,1E-1]

num_trial = 20 


def vrfunc(a,r): 
    if r >= 0 and r < 1: 
        vr = -a*(-3*r+9/4*r**2)
    elif r >= 1 and r < 2: 
        vr = 3/4*a*(2-r)**2
    else:
        vr = 0
    return vr


def compute_Ek(a,mi,Nsn): 
    Ek_parts = []
    for ipart in range(int(Nsn)):
        ri = rd.uniform(0,2) 
        Ek_parts.append(1./2.*mi*umass*(vrfunc(a,ri))**2) 
    Ek = np.sum(Ek_parts)
    return Ek 


def main(): 
    
    a_table = [] 
    Nsn_table = [] 
    mi_table = [] 
    Ek_table = [] 
    totentry = len(a_scan)*len(Nsn_scan)*len(mi_scan) 
    entry = 0
    
    for a in a_scan: 
        for Nsn in Nsn_scan: 
            for mi in mi_scan: 
                
                Ek_trials = []
                for i in range(num_trial): 
                    Ek = compute_Ek(a,mi,Nsn) 
                    Ek_trials.append(Ek) 
                meanEk = np.mean(Ek_trials)  
                
                a_table.append(a) 
                Nsn_table.append(int(np.round(Nsn,0))) 
                mi_table.append(mi*umass) 
                Ek_table.append(meanEk) 
                
                entry += 1
                if (entry/totentry*100)%1 == 0:
                    print('loading',entry/totentry*100,'%')
                
    vprofile_table = np.vstack([a_table,Nsn_table,mi_table,Ek_table]).T 
    print(pd.DataFrame(vprofile_table,columns=['a','Nsn','m','Ek']))
    np.savetxt('vprofile_sn.txt',vprofile_table,fmt='%.5e')
        
    return 



if __name__ == '__main__': 
    main()



# Phantom - Fetch from table 
# Loop through whole table to locate the entry where all vals agree 
# within tolerance = slightly bigger than param interval in table 







































