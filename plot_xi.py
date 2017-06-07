#!/usr/bin/env python2.7

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import datetime
import os


def main():
    # matplotlib options
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', family='serif')
    matplotlib.rc('lines', linewidth=3)

    # read in data
    data_python = np.genfromtxt('xi_nu.txt')
    data_fortran = np.genfromtxt('xi_nu_fortran.txt')

    # make figure
    fig = plt.figure()
    ax = plt.gca()

    nu_p = data_python[:, 0]
    xi_re_p = data_python[:, 1]
    xi_im_p = data_python[:, 2]

    nu_f = data_fortran[:, 0]
    xi_re_f = data_fortran[:, 1]
    xi_im_f = data_fortran[:, 2]
    xi_abs_f = xi_re_f**2 + xi_im_f**2

    # plot
    # ax.plot(nu_p, xi_re_p, color='g', linestyle='-', alpha=0.5, label='python')
    # ax.plot(nu_p, xi_im_p, color='g', linestyle='--', alpha=0.5)
    ax.plot(nu_f, xi_re_f, color='b', linestyle='-',
            alpha=0.5, label='real')
    ax.plot(nu_f, xi_im_f, color='b', linestyle='--', alpha=0.5, label='imag')
    # ax.plot(nu_f, xi_abs_f, color='b', linestyle='-', alpha=0.5, label='abs')

    # labels
    ax.set_xlabel(r'$\nu$ [MHz]')
    ax.set_ylabel(r'$\Xi(\nu)$')
    leg = ax.legend(loc=0)

    # save
    date = datetime.date.today().strftime('%y%m%d')
    plots_dir = ''.join(['plots/', date])
    if not os.path.isdir(plots_dir):
        os.makedirs(plots_dir)
    output = ''.join([plots_dir, '/xi_comp.pdf'])
    print("Saving {}...".format(output))
    plt.savefig(output, bbox_inches='tight')
    plt.close(fig)

    return


if __name__ == '__main__':
    main()
