import numpy as np 
import matplotlib.pyplot as plt 

# cropped voronoi sites 
x_before,y_before,z_before = np.loadtxt('before_mod_grid.txt',unpack=True)
x_after,y_after,z_after = np.loadtxt('after_mod_grid.txt',unpack=True)

pos_before = np.stack((x_before,y_before,z_before),axis=1)
pos_after = np.stack((x_after,y_after,z_after),axis=1)

print(len(pos_before))
print(len(pos_after))

minx = min(x_before)
maxx = max(x_before)
miny = min(y_before)
maxy = max(y_before)

list_before = zip(x_before,y_before,z_before)
list_after = zip(x_after,y_after,z_after)
set_before = set(list_before)
set_after = set(list_after)
sites_eaten = list(sorted(set_before - set_after))
sites_created = list(sorted(set_after - set_before))
sites_unchanged = list(sorted(set_before & set_after))

print('eaten sites',len(sites_eaten))
print('created sites',len(sites_created)) 
print('unchanged sites',len(sites_unchanged))
# Note - might not be the same as output from phantom since some groups = new site = same original particle 

nslice = 30
minz = min(z_before)
maxz = max(z_before)
dz = (maxz-minz)/nslice 

for i in range(nslice):
    lowzbound = minz + i*dz
    upzbound = minz + (i+1)*dz

    xyslice_eaten = []
    for isite in range(len(sites_eaten)):
        if (sites_eaten[isite][2] > lowzbound and sites_eaten[isite][2] < upzbound):
            xyslice_eaten.append((sites_eaten[isite][0],sites_eaten[isite][1]))

    xyslice_created = []
    for isite in range(len(sites_created)):
        if (sites_created[isite][2] > lowzbound and sites_created[isite][2] < upzbound):
            xyslice_created.append((sites_created[isite][0],sites_created[isite][1]))

    xyslice_unchanged = []
    for isite in range(len(sites_unchanged)):
        if (sites_unchanged[isite][2] > lowzbound and sites_unchanged[isite][2] < upzbound):
            xyslice_unchanged.append((sites_unchanged[isite][0],sites_unchanged[isite][1]))

    plt.figure()
    if (len(xyslice_unchanged) > 0):
        plt.scatter(np.array(xyslice_unchanged)[:,0],np.array(xyslice_unchanged)[:,1],color='black',s=0.1,marker='.')
    if (len(xyslice_eaten) > 0):
        plt.scatter(np.array(xyslice_eaten)[:,0],np.array(xyslice_eaten)[:,1],color='blue',s=0.1,marker='.')
    if (len(xyslice_created) > 0):
        plt.scatter(np.array(xyslice_created)[:,0],np.array(xyslice_created)[:,1],color='red',s=0.1,marker='.')
    plt.xlim([-0.5,2.5])
    plt.ylim([0,3])
#    plt.xlim([minx,maxx])
#    plt.ylim([miny,maxy])
    plt.savefig('slice_'+str(i)+'_voronoi_sites.png',dpi=300)










