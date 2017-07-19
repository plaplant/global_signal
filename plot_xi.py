#!/usr/bin/env python2.7

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import datetime
import os
import h5py


def main():
    # matplotlib options
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', family='serif')
    matplotlib.rc('lines', linewidth=1)

    # read in data
    f = h5py.File("Output/PAPER/beam_zenith/xi_nu_phi.hdf5", 'r')
    dgrp = f["/Data"]
    phi = np.asarray(dgrp["phi"])
    nu = np.asarray(dgrp["nu"])
    xi = np.asarray(dgrp["xi"])

    # make figure
    fig = plt.figure()
    ax = plt.gca()

    # add data
    add_data_to_plot(ax, nu, phi, xi)

    # labels and things
    ax.set_xlabel(r'$\nu$ [MHz]')
    ax.set_ylabel(r'$\Xi(\nu)$')
    ax.set_xlim((100, 200))
    #leg = ax.legend(loc=0)

    # save
    date = datetime.date.today().strftime('%y%m%d')
    plots_dir = ''.join(['plots/', date])
    if not os.path.isdir(plots_dir):
        os.makedirs(plots_dir)
    output = ''.join([plots_dir, '/paper_xi_comp.pdf'])
    print("Saving {}...".format(output))
    plt.savefig(output, bbox_inches='tight')
    plt.close(fig)

    return


def add_data_to_plot(ax, nu, phi, data):
    # make line color cycle thorugh viridis colormap
    nphi = len(phi)
    color_idx = np.linspace(0, 1, nphi)
    for i, c in enumerate(color_idx):
        x = nu
        yre = np.real(data[i, :])
        yim = np.imag(data[i, :])
        ax.plot(x, yre, color=plt.cm.viridis(c), linestyle='-', alpha=0.5)
        ax.plot(x, yim, color=plt.cm.viridis(c), linestyle='--', alpha=0.5)

    return


if __name__ == '__main__':
    main()
