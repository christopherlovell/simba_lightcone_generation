import sys

import h5py
import numpy as np
from astropy.cosmology import z_at_value
import astropy.units as u

from simba import simba


sb = simba()
verbose = True

lc_fname = 'output/lightcone_%03d.h5'

N_lcs = int(sys.argv[1])
area =  float(sys.argv[2]) # degree ^ 2
z_min = float(sys.argv[3])
z_max = float(sys.argv[4])

if N_lcs < 1: raise ValueError('Need at least one lightcone to generate!')

outs = np.loadtxt(sb.output_file)
snaps = sb.lightcone_snaps
zeds = np.array([1./outs[int(snap)] - 1 for snap in snaps])

## ---- check area can be covered by box
cs = sb.get_caesar(snaps[0])
L = cs.simulation.boxsize.to('Mpccm').value

zeds_mask = (zeds > z_min) & (zeds < z_max)
L_unit = sb.cosmo.kpc_comoving_per_arcmin(zeds[zeds_mask]).to('Mpc / degree')
A = (L_unit * area**0.5).value

if np.any(A > L):
    raise ValueError('Specified area too large for simulation box (lateral tiling not yet implemented)')

## ---- start lightcone creation
lc_out = {str(_lc): {} for _lc in np.arange(N_lcs)}

_N_all = 0
for i,snap in enumerate(snaps[zeds_mask]):
    z = 1./outs[int(snap)] - 1
    print("\nz:",z,snap)

    cs = sb.get_caesar(snap)
    a = cs.simulation.scale_factor    

    L_unit = sb.cosmo.kpc_comoving_per_arcmin(z).to('Mpc / degree')
    A = (L_unit * area**0.5).value
    z_offset = sb.cosmo.comoving_distance(z).value - L/2
    
    coods_cMpc = np.array([g.pos.to('Mpccm').value for g in cs.galaxies])

    _N_all += len(coods_cMpc)
    
    if verbose: print("Generating lightcone selection...")
    lc_mask = np.ones(len(coods_cMpc),dtype=bool)

    for _lc in np.arange(N_lcs):
        lc_out[str(_lc)][snap] = {}
        if verbose: print("N lightcone:", _lc)

        i = np.random.randint(0,3)  # randomly choose axes
        j = i
        while j == i: j = np.random.randint(0,3)
        k = np.where((np.arange(0,3) != i) & (np.arange(0,3) != j))[0][0]
    
        xmin,ymin = np.random.rand(2) * (L-A)   # get minimum box mask coordinate
        if verbose: print("xmin:",xmin, "\nymin:",ymin, "\nA:",A, "\nL:",L)
   
        lc_idx_arr = lc_mask & ((coods_cMpc[:,i] > xmin) &\
                                (coods_cMpc[:,i] < xmin+A) &\
                                (coods_cMpc[:,j] > ymin) &\
                                (coods_cMpc[:,j] < ymin+A))
        
        if verbose: print("N(lightcone cut):",np.sum(lc_idx_arr))
    
        lc_idx_arr = np.where(lc_idx_arr)[0]
        
        _coods = coods_cMpc[lc_idx_arr]
        _coods[:,i] -= xmin
        _coods[:,j] -= ymin
        _coods = _coods[:,[i,j,k]]

        _ra = _coods[:,0] / L_unit # deg
        _dec = _coods[:,0] / L_unit # deg
        _redshift = [z_at_value(sb.cosmo.comoving_distance, _c * u.Mpc) \
                          for _c in (_coods[:,2] + z_offset)]

        lc_out[str(_lc)][snap]['index'] = lc_idx_arr # .tolist()
        lc_out[str(_lc)][snap]['RA'] = _ra.value # .tolist()
        lc_out[str(_lc)][snap]['DEC'] = _dec.value # .tolist()
        lc_out[str(_lc)][snap]['z'] = np.array(_redshift)


print("_N_all:",_N_all)

print("Writing lightcone selection to %s"%lc_fname)
for _lc in np.arange(N_lcs):
    sb.save_dict_to_hdf5(lc_out[str(_lc)], lc_fname%_lc)

