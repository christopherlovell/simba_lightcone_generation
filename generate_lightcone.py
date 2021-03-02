import sys
import argparse

import h5py
import numpy as np
from astropy.cosmology import z_at_value
import astropy.units as u

from simba import simba


sb = simba()
verbose = True

parser = argparse.ArgumentParser(description='Generate a Simba lightcone.', add_help=True,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('area', type=float, help="Survey area (deg^2)")
parser.add_argument('z_min', type=float, help="Minimum redshift")
parser.add_argument('z_max', type=float, help="Maximum redshift")

parser.add_argument('--N', type=int, default=1, action='store', 
                    help="Number of lightcones.")

parser.add_argument('-f', '--filename', 
                    action='store', 
                    dest='filename',
                    type=str, 
                    #required=False, 
                    default='output/halo_lightcone_%03d.h5', 
                    help="Output filename (Full path from current directory. \
                          Must include string formatting to iterate over lightcone number)")

parser.add_argument("--overwrite", help="Overwrite existing datasets", 
                    action='store_true')

parser.add_argument("--halos", help="Output halos, rather than galaxies.", 
                    action='store_true')

parser.add_argument("-v", "--verbose", help="Increase output verbosity", 
                    action='store_true')

args = parser.parse_args()

verbose = args.verbose
area = args.area
N_lcs = args.N
z_min = args.z_min
z_max = args.z_max
lc_fname = args.filename
overwrite = args.overwrite
halo_flag = args.halos
lc_fname = 'output/halo_lightcone_%03d.h5'

# N_lcs = int(sys.argv[1])
# area =  float(sys.argv[2]) # degree ^ 2
# z_min = float(sys.argv[3])
# z_max = float(sys.argv[4])

if N_lcs < 1: raise ValueError('Need at least one lightcone to generate!')

outs = np.loadtxt(sb.output_file)
snaps = sb.lightcone_snaps
zeds = np.array([1./outs[int(snap)] - 1 for snap in snaps])

## ---- check area can be covered by box
cs = sb.get_caesar(snaps[0])
L = cs.simulation.boxsize.to('Mpccm').value

zeds_mask = (zeds >= z_min) & (zeds <= z_max)
L_unit = sb.cosmo.kpc_comoving_per_arcmin(zeds[zeds_mask]).to('Mpc / degree')
A = (L_unit * area**0.5).value

if np.any(A > L):
    raise ValueError('Specified area too large for simulation box (lateral tiling not yet implemented)')

## ---- start lightcone creation
lc_out = {str(_lc): {} for _lc in np.arange(N_lcs)}

_N_all = 0
# z_prev = 0.
A_A = 0.
for i,(snapA,snapB) in enumerate(zip(snaps[zeds_mask],snaps[zeds_mask][1:])):
    
    z   = 1./outs[int(snapA)] - 1
    z_B = 1./outs[int(snapB)] - 1
    print("\nz:",z,snapA)

    cs = sb.get_caesar(snapA)
    a = cs.simulation.scale_factor    

    L_unit = sb.cosmo.kpc_comoving_per_arcmin(z_B).to('Mpc / degree')
    A = (L_unit * area**0.5).value

    z_offset = sb.cosmo.comoving_distance(z).value
    
    if halo_flag:
        coods_cMpc = np.array([h.pos.to('Mpccm').value for h in cs.halos])
    else:
        coods_cMpc = np.array([g.pos.to('Mpccm').value for g in cs.galaxies])

    _N_all += len(coods_cMpc)
    
    if verbose: print("Generating lightcone selection...")
    # lc_mask = np.ones(len(coods_cMpc),dtype=bool)

    for _lc in np.arange(N_lcs):
        lc_out[str(_lc)][snapA] = {}
        if verbose: print("N lightcone:", _lc)

        i = np.random.randint(0,3)  # randomly choose axes
        j = i
        while j == i: j = np.random.randint(0,3)
        k = np.where((np.arange(0,3) != i) & (np.arange(0,3) != j))[0][0]
    
        xmin,ymin = np.random.rand(2) * (L-A)   # get minimum box mask coordinate
        if verbose: print("xmin:",xmin, "\nymin:",ymin, "\nA:",A, "\nL:",L)
  
        # lightcone 'frustum' angle
        theta = np.arctan((A - A_A) / (2*L))
        # dx == distance along cood x between top and bottom of fustrum edge
        dx = np.abs(L - coods_cMpc[:,k]) * np.tan(theta)

        lc_idx_arr = ((coods_cMpc[:,i] > (xmin + dx)) &\
                      (coods_cMpc[:,i] < ((xmin+A) - dx)) &\
                      (coods_cMpc[:,j] > (ymin + dx)) &\
                      (coods_cMpc[:,j] < ((ymin+A) - dx)))
        
        if verbose: print("N(lightcone cut):",np.sum(lc_idx_arr))
    
        lc_idx_arr = np.where(lc_idx_arr)[0]        
        _coods = coods_cMpc[lc_idx_arr]

        _frac = np.abs(_coods[:,i] - xmin - (A/2)) / ((A/2) - dx[lc_idx_arr])
        _coods[:,i] = _frac * ((A * u.Mpc) / L_unit)
        
        _frac = np.abs(_coods[:,j] - ymin - (A/2)) / ((A/2) - dx[lc_idx_arr])
        _coods[:,j] = _frac * ((A * u.Mpc) / L_unit)

        _coods = _coods[:,[i,j,k]]

        _ra  = _coods[:,0]
        _dec = _coods[:,1]
        _redshift = np.array([z_at_value(sb.cosmo.comoving_distance, _c * u.Mpc) \
                              for _c in (_coods[:,2] + z_offset)])

        # lc_out[str(_lc)][snapA]['index'] = lc_idx_arr
        # lc_out[str(_lc)][snapA]['RA'] = _ra 
        # lc_out[str(_lc)][snapA]['DEC'] = _dec 
        # lc_out[str(_lc)][snapA]['z'] = np.array(_redshift)

        with h5py.File(lc_fname%_lc,'a') as h5file:
            h5file.require_group(snapA)

        sb.create_dataset(lc_fname%_lc, lc_idx_arr, 'index', group=snapA, overwrite=overwrite)
        sb.create_dataset(lc_fname%_lc, _ra, 'RA', group=snapA, overwrite=overwrite)
        sb.create_dataset(lc_fname%_lc, _dec, 'DEC', group=snapA, overwrite=overwrite)
        sb.create_dataset(lc_fname%_lc, _redshift, 'z', group=snapA, overwrite=overwrite)


    A_A = A


print("_N_all:",_N_all)

# print("Writing lightcone selection to %s"%lc_fname)
# for _lc in np.arange(N_lcs):
#     sb.save_dict_to_hdf5(lc_out[str(_lc)], lc_fname%_lc)

