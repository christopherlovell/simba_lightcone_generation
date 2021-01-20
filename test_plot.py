import json

import matplotlib.pyplot as plt


lc_fname = 'output/lightcone_%03d.json'
_lc = 0

with open(lc_fname%_lc, 'r') as fp:
    lc_out = json.load(fp)


for snap in lc_out.keys():     
    zeds = lc_out[str(snap)]['z'] 
    RA = lc_out[str(snap)]['RA'] 
    DEC = lc_out[str(snap)]['DEC'] 
    plt.scatter(zeds,RA) 

plt.xlabel('z')
plt.ylabel('R.A. (deg)')
plt.show()   

