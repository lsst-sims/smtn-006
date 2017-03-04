"""
This script reads in the HEALPIX maps generated by

stellar_density_get_arrays.py
stellar_density_control_arrays.py

and creates plots of them.
"""

from __future__ import with_statement

import matplotlib
matplotlib.use('Agg')

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

import time

import os

import argparse


# user-specified resolution.  Just make sure that it matches with
# stellar_density_get_arrays.py and stellar_density_control_arrays.py
NSIDE = 64


def plot_density(input_dir_list, mag, i_bound, i_fig, rows=3, cols=2, title_list=None):

    n_pix = hp.nside2npix(NSIDE)

    margins = (0.01, 0.01, 0.04, 0.01)

    dtype = np.dtype([('ct', int)])
    tag_list = []

    legend = None

    tag = '%s_%d' % (mag, i_bound)

    ct_arr = {}

    for input_dir in input_dir_list:
        list_of_files = os.listdir(input_dir)
        ct_arr[input_dir] = np.zeros(n_pix)
        for file_name in list_of_files:

            if not tag in file_name:
                continue

            full_name = os.path.join(input_dir, file_name)

            with open(full_name, 'r') as input_file:
                if legend is None:
                    legend = input_file.readlines()[0]
                    print legend
                else:
                    assert legend == input_file.readlines()[0]

            data = np.genfromtxt(full_name, dtype=dtype)
            assert len(data['ct']) == n_pix
            ct_arr[input_dir] += data['ct']

    c_min = None
    c_max = None
    for input_dir in input_dir_list:
        lmin = ct_arr[input_dir].min()
        lmax = ct_arr[input_dir].max()
        if c_min is None or lmin<c_min:
            c_min = lmin
        if c_max is None or lmax>c_max:
            c_max = lmax

    for input_dir, title in zip(input_dir_list, title_list):
        if title_list is None:
            fig_title = legend.strip().replace('# ','')
        else:
            fig_title = legend.strip().replace('# ','') + ' ' + title

        hp.mollview(ct_arr[input_dir], title=fig_title,cbar=False, margins=margins,
                    sub=(rows, cols, i_fig))
        hp.graticule(dpar=10, dmer=20, verbose=False)

        i_fig += 1

        ax = plt.gca()
        im = ax.get_images()[0]

        cticks = np.arange(c_min, c_max, 0.2*(c_max-c_min))
        clabels = ['%.2e' % cc for cc in cticks]

        cb = plt.colorbar(im)
        cb.set_ticks(cticks)
        cb.set_ticklabels(clabels)
        cb.set_clim(vmin=c_min, vmax=c_max)
        cb.ax.tick_params(labelsize=10)
        cb.draw_all()

    return i_fig

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, default=None)
    parser.add_argument("--output_dir", type=str, default=None)

    args = parser.parse_args()
    if args.input_dir is None:
        raise RuntimeError("Must supply input_dir")
    if args.output_dir is None:
        raise RuntimeError("Must supply output_dir")

    if os.path.exists(args.output_dir) and not os.path.isdir(args.output_dir):
        raise RuntimeError('%s is not a dir' % args.output_dir)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    fit_map_dir = os.path.join(args.input_dir, 'fit_healpix_maps')
    input_map_dir = os.path.join(args.input_dir, 'input_healpix_maps')

    with open(os.path.join(fit_map_dir, 'number_of_bounds.txt'), 'r') as in_file:
        n_bounds_fit = int(in_file.readlines()[0])

    with open(os.path.join(input_map_dir, 'number_of_bounds.txt'), 'r') as in_file:
        n_bounds_input = int(in_file.readlines()[0])

    assert n_bounds_fit == n_bounds_input
    print 'n_bounds ',n_bounds_fit

    mag_list = ['u', 'g', 'r', 'i', 'z', 'y']

    rows = 3
    cols = 2

    for mag in mag_list:

        plt.figsize = (30,30)
        i_plot=0
        i_fig = 1
        fig_name = os.path.join(args.output_dir, '%s_fig_%d.png' % (mag, i_plot))

        for ib in range(n_bounds_fit):

            i_fig = plot_density([input_map_dir, fit_map_dir], mag, ib, i_fig,
                                 title_list=['(input)', '(fit)'], rows=rows, cols=cols)

            if i_fig>rows*cols:
                plt.savefig(fig_name)
                plt.close()
                i_plot += 1
                fig_name = os.path.join(args.output_dir, '%s_fig_%d.png' % (mag, i_plot))
                i_fig = 1
                plt.figsize=(30, 30)

        if i_fig>1:
            plt.savefig(fig_name)
            plt.close()
            i_plot += 1
            i_fig = 1
