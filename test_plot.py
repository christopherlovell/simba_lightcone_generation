import json
import matplotlib.pyplot as plt
from simba import simba

sb = simba()

lc_fname = 'output/test_lightcone_%03d.h5'
_lc = 0

lc_out = sb.load_dict_from_hdf5(lc_fname%_lc)

fig, (ax1,ax2) = plt.subplots(2,1,figsize=(5,8))

for snap in lc_out.keys():     
    zeds = lc_out[str(snap)]['z'] 
    RA = lc_out[str(snap)]['RA'] 
    DEC = lc_out[str(snap)]['DEC'] 
    ax1.scatter(zeds,RA) 
    ax2.scatter(zeds,DEC) 

ax2.set_xlabel('z')
ax1.set_ylabel('R.A. (deg)')
ax2.set_ylabel('Dec. (deg)')

plt.show()   

