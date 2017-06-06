#!/usr/bin/env python2.7

from __future__ import print_function, division, absolute_import
import numpy as np
import healpy as hp


def test_raw_data():
    datafile = '/data4/paper/zionos/polskysim/IonRIME/jones_save/HERA_NicCST/ijonesband_100-200mhz_nfreq201_nside128.npy'
    data = np.load(datafile, mmap_mode='r')

    # check certain indices
    i = 0
    j = 0
    k = 147714
    l = 150
    z = data[l, k, j, i]
    print("input_maps[{:d},{:d},{:d},{:d}]: {:16.8e} {:16.8e}j".format(
        l, k, j, i, z.real, z.imag))

    i = 0
    j = 1
    k = 147714
    l = 125
    z = data[l, k, j, i]
    print("input_maps[{:d},{:d},{:d},{:d}]: {:16.8e} {:16.8e}j".format(
        l, k, j, i, z.real, z.imag))

    return


def test_ii_map():
    datafile = '/data4/paper/zionos/polskysim/IonRIME/jones_save/HERA_NicCST/ijonesband_100-200mhz_nfreq201_nside128.npy'
    data = np.load(datafile, mmap_mode='r')

    # compute specific pixels
    ip = 147712
    inu = 150
    p = 0
    for i in range(2):
        for j in range(2):
            p += np.abs(data[inu, ip, i, j])**2
    print("ii_maps[{:d},{:d}]: {:16.8e}".format(inu, ip, p))

    ip = 147677
    inu = 100
    p = 0
    for i in range(2):
        for j in range(2):
            p += np.abs(data[inu, ip, i, j])**2
    print("ii_maps[{:d},{:d}]: {:16.8e}".format(inu, ip, p))

    return


def test_alm():
    datafile = '/data4/paper/zionos/polskysim/IonRIME/jones_save/HERA_NicCST/ijonesband_100-200mhz_nfreq201_nside128.npy'
    data = np.load(datafile, mmap_mode='r')
    NSIDE = 128
    LMAX = 3 * NSIDE - 1

    # compute map for whole frequency
    inu = 97
    mp = np.zeros((data.shape[1]))
    for i in range(2):
        for j in range(2):
            mp[:] += np.abs(data[inu, :, i, j])**2

    # compare with map from fortran
    mp_f = np.genfromtxt('test_map.txt')
    print("python and fortran maps agree: ", np.allclose(mp_f, mp))

    # take alm transform
    alm = hp.sphtfunc.map2alm(mp, lmax=LMAX)

    # print a few specific values
    l = 0
    m = 0
    idx = hp.sphtfunc.Alm.getidx(LMAX, l, m)
    a = alm[idx]
    print("a_{{{:d},{:d}}}: {:16.8e} {:16.8e}j".format(l, m, a.real, a.imag))
    l = 7
    m = 3
    idx = hp.sphtfunc.Alm.getidx(LMAX, l, m)
    a = alm[idx]
    print("a_{{{:d},{:d}}}: {:16.8e} {:16.8e}j".format(l, m, a.real, a.imag))
    l = 128
    m = 99
    idx = hp.sphtfunc.Alm.getidx(LMAX, l, m)
    a = alm[idx]
    print("a_{{{:d},{:d}}}: {:16.8e} {:16.8e}j".format(l, m, a.real, a.imag))

    return


if __name__ == '__main__':
    test_raw_data()
    test_ii_map()
    test_alm()
