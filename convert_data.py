#!/usr/bin/env python2.7

from __future__ import print_function, division, absolute_import
import numpy as np
import h5py


def convert_data(infile, outfile):
    # load data
    data = np.load(infile)

    # make hdf5 file
    f = h5py.File(outfile, "w")
    f.attrs['Nfreq'] = 201
    f.attrs['Nside'] = 128
    f.attrs['band_low'] = 100
    f.attrs['band_high'] = 200
    f.attrs['band_unit'] = "MHz"

    # write out data
    grp = f.create_group("Data")
    ds = grp.create_dataset("MuellerMatrices", data=data)
    f.close()


if __name__ == '__main__':
    infile = '/data4/paper/zionos/polskysim/IonRIME/jones_save/HERA_NicCST/ijonesband_100-200mhz_nfreq201_nside128.npy'
    outfile = '/data4/paper/plaplant/beams/HERA_ijones.hdf5'
    convert_data(infile, outfile)
