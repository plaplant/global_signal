#!/usr/bin/env python2.7

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import healpy as hp
from numba import jit
from scipy.special import spherical_jn, sph_harm
import datetime
import os

datafile = '/data4/paper/zionos/polskysim/IonRIME/jones_save/HERA_NicCST/ijonesband_100-200mhz_nfreq201_nside128.npy'
outfile = '/home/plaplant/global_signal/xi_nu.txt'

# constants for calculation
NFREQ = 201
NSIDE = 128
NPIX = 12 * NSIDE**2
TAUH = 50e-9
LMAX = 3 * NSIDE - 1


@jit
def compute_II_map(arr):
    nf, npix, ni, nj = arr.shape
    output = np.zeros((NFREQ, NPIX))
    for i in range(ni):
        for j in range(nj):
            output[:, :] += np.abs(arr[:, :, i, j])**2

    return output


def calc_xi_nu(alm, nu, tauh, lmax):
    """
    Calculate xi(nu) for a particular frequency given spherical harmonic coefficients

    Args:
       alm (array) -- array of HEALPix spherical harmonic transform coefficients.
       nu (float) -- the frequency being calculated. Units of s^{-1}.
       tauh (float) -- the speed-of-light delay of the target baseline. Units of s.
       lmax (int) -- the maximum l-value of the HEALPix transform.

    Output:
       xi (complex) -- xi(nu) for the given harmonic coefficients
    """
    # initialize
    xi = 0j

    # loop over l and m values
    for il in range(lmax + 1):
        arg = 2 * np.pi * tauh * nu
        jl = spherical_jn(il, arg)

        for im in range(il + 1):
            idx = hp.sphtfunc.Alm.getidx(lmax, il, im)
            prefac = (1j)**(3 * il + 2 * im)
            ylm = sph_harm(im, il, 0., 0.)
            a = alm[idx]

            # add to running total
            if np.isnan(prefac) or np.isnan(a) or np.isnan(jl) or np.isnan(ylm):
                continue
            xi += prefac * a * jl * ylm
            if im > 0:
                # for negative m, we have:
                #   a_{l, -m} = a_{l, m}^*
                #   Y_l^{-m} = (-1)**m * (Y_l^m)^*
                prefac *= (-1)**im
                if np.isnan(prefac) or np.isnan(np.conj(a)) or np.isnan(np.conj(ylm)):
                    continue
                xi += prefac * np.conj(a) * jl * np.conj(ylm)

    return np.complex_(xi)


def make_II_map(data):
    """
    Make a healpy map of the I -> I' beam component at all frequencies.

    Args:
       data (str) -- the location of the raw data file containing the Jones matrices
       freq (float) -- the frequency (in MHz) of the desired output

    Output:
       mp (numpy array, size(Nfreq, Npix)) -- an nside=128 HealPIX map containing the I -> I' beam
    """

    ijones = np.load(data, mmap_mode='r')

    mp = compute_II_map(ijones)

    return mp


def compute_xi(mp):
    """
    Calculate the xi(nu) function for the global signal paper.

    Args:
       mp (numpy array) -- a HealPIX map containing the beam. Assumes the map
          contains 201 frequencies evenly spaced between 100 and 200 MHz

    Output:
       xi -- an array as a function of frequency of response of an interferometer
          to the global signal
    """
    # initialize frequency array
    nu = np.linspace(100, 200, num=NFREQ, endpoint=True)

    # convert from MHz -> Hz
    nu = nu * 1e6

    # initialize output array
    xi_nu = np.zeros(NFREQ, dtype=np.complex_)

    # loop over frequencies
    for i in range(NFREQ):
        print("iteration {:d}".format(i))
        # get healpix map and convert to a_lm's
        mp1 = mp[i, :]
        alm = hp.sphtfunc.map2alm(mp1, lmax=LMAX)
        print("spherical transform done")

        # iterate over l,m values
        xi_nu[i] = calc_xi_nu(alm, nu[i], TAUH, LMAX)

    # add overall normalization
    xi_nu *= np.sqrt(4 * np.pi)

    return nu, xi_nu


def main():
    try:
        # read in the pre-computed data
        data = np.genfromtxt(outfile)
        nu = data[:, 0]
        xi_nu = data[:, 1] + 1j * data[:, 2]
    except IOError:
        # compute the Mueller I -> I' map
        print("Computing I -> I' map...")
        mp = make_II_map(datafile)

        # compute xi(nu) for the map
        print("Computing xi(nu)...")
        nu, xi_nu = compute_xi(mp)

        # save the data
        print("Saving data...")
        output = np.empty((len(nu), 3))
        output[:, 0] = nu
        output[:, 1] = np.real(xi_nu)
        output[:, 2] = np.imag(xi_nu)
        np.savetxt(outfile, output, header=' nu    xi(nu)')

    # plotting options
    #matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', family='serif')
    matplotlib.rc('lines', linewidth=3)

    # plot it
    fig = plt.figure()
    ax = plt.gca()

    # convert to MHz
    nu /= 1e6

    # take absolute value
    xi_nu = np.abs(xi_nu)

    # plot
    ax.plot(nu, xi_nu, linestyle='-', color='r')

    # labels
    ax.set_xlabel('nu [MHz]')
    ax.set_ylabel('Xi(nu)')

    # save plot
    date = datetime.date.today().strftime('%y%m%d')
    plots_dir = ''.join(['/home/plaplant/global_signal/plots/', date])
    if not os.path.isdir(plots_dir):
        os.makedirs(plots_dir)
    output = ''.join([plots_dir, '/xi_nu.pdf'])
    print("Saving {}...".format(output))
    plt.savefig(output, bbox_inches='tight')
    plt.close(fig)

    return


if __name__ == '__main__':
    main()
