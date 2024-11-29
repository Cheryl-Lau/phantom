
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

vel, temp = np.loadtxt('vel_temp_cloud_20_10_clsink139_v2_06456.dat',unpack=True)

plt.scatter(temp,vel,s=1,color='black')

plt.xlabel('T [K]')
plt.ylabel('v [cm/s]')
plt.xscale('log')
plt.yscale('log')

plt.show()
