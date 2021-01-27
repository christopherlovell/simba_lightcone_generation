# Simba lightcone generation

## Requirements
- numpy
- astropy
- caesar

## Running

`python generate_lightcone.py N A zmin zmax`

where `N` is the number of lightcones you wish to generate, `A` is the sky area in deg<sup>2</sup> (must be small enough that it covers the whole box, as lateral tiling is not yet implemented; 0.5 deg<sup>2</sup> works), and `zmin` and `zmax` are the lower and upper redshift limits.

## Output

Each lightcone is saved in its own `json` file in the `output` folder. These are arranged by a top folder key for each snapshot, within which are lists for RA, DEC, index and redshift ('z'). The index refers back to the caesar galaxy lists.
