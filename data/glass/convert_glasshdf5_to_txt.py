
import numpy as np
import pandas as pd
import h5py

f = h5py.File('glass256.hdf5',"r")


#for key in f["/PartType0/Coordinates"]:
#   print(key)

coords = f.get("/PartType0/Coordinates")
coords = np.array(coords)

meanx = np.mean(coords[:,0])
meany = np.mean(coords[:,1])
meanz = np.mean(coords[:,2])

print(meanx,meany,meanz)

coords[:,0] = coords[:,0] - meanx
coords[:,1] = coords[:,1] - meany
coords[:,2] = coords[:,2] - meanz

meanx = np.mean(coords[:,0])
meany = np.mean(coords[:,1])
meanz = np.mean(coords[:,2])

print(meanx,meany,meanz)



h = f.get("/PartType0/SmoothingLength")
h = np.array(h)

df = pd.DataFrame(coords, columns =['x', 'y','z']) 
df['h'] = h

glass = df.to_numpy()

print(len(h))

np.savetxt('glassCube_256.dat',glass)




