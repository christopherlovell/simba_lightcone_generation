import numpy as np
import caesar

from astropy.cosmology import Planck15 as cosmo
from astropy import units as u


class simba:
    def __init__(self):

        self.lightcone_snaps = np.array([str(s).zfill(3) for s in np.arange(20,145,2)[::-1]])

        self.sim_directory='/cosma7/data/dp104/dc-dave2/sim/m100n1024/s50j7k/'
        # self.cs_directory=self.sim_directory+'Groups/'#'Groups_old/caesar_old/'
        self.cs_directory=self.sim_directory+'Groups_old/caesar_old/'
        self.output_file='/cosma7/data/dp104/dc-dave2/sim/m100n1024/s50/outputs_boxspace50.txt'
        self.cosmo = cosmo

    def get_sim_file(self,snap,snap_str="snap_m100n1024_%s.hdf5"):
        #return self.sim_directory+('snap_m100n1024_%s.hdf5'%snap)
        return self.sim_directory+(snap_str%snap)

    def get_caesar(self,snap,fname='m100n1024_%s.hdf5',verbose=False):
        fname = self.cs_directory+(fname%snap)
        # return  caesar.quick_load(fname)
        return  caesar.load(fname)

