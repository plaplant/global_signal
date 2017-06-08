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
    data1 = np.genfromtxt('Output/beam_default/b_zenith/xi_nu.txt')
    data2 = np.genfromtxt('Output/beam_default/b_horizon/xi_nu.txt')
    data3 = np.genfromtxt('Output/beam_default/b_horizon2/xi_nu.txt')
    data4 = np.genfromtxt('Output/beam_default/b_horizon3/xi_nu.txt')
    data5 = np.genfromtxt('Output/beam_zenith/b_zenith/xi_nu.txt')
    data6 = np.genfromtxt('Output/beam_zenith/b_horizon/xi_nu.txt')
    data7 = np.genfromtxt('Output/beam_zenith/b_horizon2/xi_nu.txt')
    data8 = np.genfromtxt('Output/beam_zenith/b_horizon3/xi_nu.txt')

    # make figure
    fig = plt.figure()
    ax = plt.gca()

    #add_data_to_plot(ax, data1, color='g', label='default zenith')
    #add_data_to_plot(ax, data2, color='b', label='default horizon')
    #add_data_to_plot(ax, data3, color='y', label='default horizon2')
    #add_data_to_plot(ax, data4, color='g', label='default horizon3')
    #add_data_to_plot(ax, data5, color='c', label='zenith zenith')
    add_data_to_plot(ax, data6, color='m', label='zenith horizon')
    add_data_to_plot(ax, data7, color='y', label='zenith horizon2')
    add_data_to_plot(ax, data8, color='k', label='zenith horizon3')

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


def add_data_to_plot(ax, data, color='k', label=None):
    x = data[:, 0]
    re = data[:, 1]
    im = data[:, 2]

    ax.plot(x, re, color=color, linestyle='-', alpha=0.5, label=label)
    ax.plot(x, im, color=color, linestyle='--', alpha=0.5)

    return


if __name__ == '__main__':
    main()
