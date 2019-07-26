#!/usr/bin/env python

#from ROOT import *
import uproot
import h5py
import numpy as np
import sys, os

#fin = TFile(ifName)
#tree = fin.Get("Events")
fin = uproot.open(sys.argv[1])
tree = fin['Events']

brNames = ["run", "event", "weight",
           "MET_pt",
           "hTrck_pt", "hEcal_pt", "hHcal_pt",
           "hTrck_n", "hEcal_n", "hHcal_n"]
#pd = tree.pandas.df(["run", "event", "weight",
#                     "MET_pt",
#                     "hTrck_pt", "hEcal_pt", "hHcal_pt",
#                     "hTrck_n", "hEcal_n", "hHcal_n"])
#pd.to_hdf(sys.argv[2], 's', complevel=9, mode='w', append=False, complib='bzip2')

with h5py.File(sys.argv[2], 'w') as fout:
    dt = h5py.special_dtype(vlen=np.float32)
    tout = fout.create_group("Events")
    for brName in brNames:
        data = np.array(tree[brName].array())
        print data.shape
        if type(data[0]) == np.ndarray:
            tout.create_dataset(brName, data.shape, data=data, chunks=True, 
                                compression='gzip', compression_opts=9)
        else:
            tout.create_dataset(brName, data=data, chunks=True)
