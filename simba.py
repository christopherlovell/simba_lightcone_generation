import numpy as np
import h5py

import caesar

from astropy.cosmology import Planck15 as cosmo
from astropy import units as u


class simba:
    def __init__(self):

        self.lightcone_snaps = np.array([str(s).zfill(3) for s in np.arange(1,152,2)[::-1]])

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

    def _check_hdf5(self, fname, obj_str):
        with h5py.File(fname, 'a') as h5file:
            if obj_str not in h5file:
                return False
            else:
                return True

    def create_dataset(self, fname, values, name, group='/', overwrite=False,
                       dtype=np.float64, desc = None, unit = None, verbose=False):

        shape = np.shape(values)

        if self._check_hdf5(fname, group) is False:
            raise ValueError("Group does not exist")
            return False

        try:
            with h5py.File(fname, mode='a') as h5f:

                if overwrite:
                    if verbose: print('Overwriting data in %s/%s'%(group,name))
                    if self._check_hdf5(fname, "%s/%s"%(group,name)) is True:
                        grp = h5f[group]
                        del grp[name]


                dset = h5f.create_dataset("%s/%s"%(group,name), shape=shape,
                                           maxshape=(None,) + shape[1:],
                                           dtype=dtype, 
                                           #compression=self.compression,
                                           data=values)

                if desc is not None:
                    dset.attrs['Description'] = desc
                if unit is not None:
                    dset.attrs['Units'] = unit

        except Exception as e:
            print("Oh! something went wrong while creating {}/{} or it already exists.\
                   \nNo value was written into the dataset.".format(group, name))
            print (e)
            # sys.exit
    

    def save_dict_to_hdf5(self, dic, filename):
        """
        ....
        """
        with h5py.File(filename, 'a') as h5file:
            self.recursively_save_dict_contents_to_group(h5file, '/', dic)
    
    def recursively_save_dict_contents_to_group(self, h5file, path, dic):
        """
        ....
        """
        for key, item in dic.items():
            if isinstance(item, (np.ndarray, np.int64, np.float64, str, bytes)):
                h5file[path + key] = item
            elif isinstance(item, dict):
                self.recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
            else:
                raise ValueError('Cannot save %s type'%type(item))


    def load_dict_from_hdf5(self, filename):
        """
        ....
        """
        with h5py.File(filename, 'r') as h5file:
            return self.recursively_load_dict_contents_from_group(h5file, '/')
    
    def recursively_load_dict_contents_from_group(self, h5file, path):
        """
        ....
        """
        ans = {}
        for key, item in h5file[path].items():
            if isinstance(item, h5py._hl.dataset.Dataset):
                ans[key] = item.value
            elif isinstance(item, h5py._hl.group.Group):
                ans[key] = self.recursively_load_dict_contents_from_group(h5file, path + key + '/')
        return ans 
