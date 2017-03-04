from __future__ import with_statement
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import argparse
import os
import numpy as np

def populate_grid(input_dir, mag1, grid):
    list_of_files = os.listdir(input_dir)

    dtype = np.dtype([('x', str, 30), ('y', str, 30), ('ct', int)])

    for file_name in list_of_files:
        if file_name.startswith(mag1):
            data = np.genfromtxt(os.path.join(input_dir, file_name), dtype=dtype)
            for xx, yy, cc in zip(data['x'], data['y'], data['ct']):
                if xx not in grid:
                    grid[xx] = {}
                if yy not in grid[xx]:
                    grid[xx][yy] = cc
                else:
                    grid[xx][yy] += cc

    return None

def plot_grid(grid, xlabel, ylabel, title, i_fig,
              xbounds=None, ybounds=None, rows=1, cols=2, cbounds=None):
    x_arr = []
    y_arr = []
    ct_arr = []
    for xx in grid:
        for yy in grid[xx]:
            x_arr.append(float(xx))
            y_arr.append(float(yy))
            ct_arr.append(grid[xx][yy])
    x_arr = np.array(x_arr)
    y_arr = np.array(y_arr)
    ct_arr = np.array(ct_arr)

    sorted_dex = np.argsort(ct_arr)

    x_arr = x_arr[sorted_dex]
    y_arr = y_arr[sorted_dex]
    ct_arr = ct_arr[sorted_dex]

    cum_arr = np.cumsum(ct_arr)

    plt.subplot(rows, cols, i_fig)
    plt.scatter(x_arr, y_arr, c=ct_arr.astype(float),
                edgecolor='',s=5,
                cmap=plt.cm.gist_ncar)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)

    if xbounds is not None:
        plt.xlim(xbounds)
    else:
        xbounds = (x_arr.min()-0.5, x_arr.max()+0.5)
        plt.xlim(xbounds)

    if ybounds is not None:
        plt.ylim(ybounds)
    else:
        ybounds = (y_arr.min()-0.5, y_arr.max()+0.5)
        plt.ylim(ybounds)

    cb = plt.colorbar()
    if cbounds is not None:
        cb.set_clim(vmin=cbounds[0], vmax=cbounds[1])

    return xbounds, ybounds, (ct_arr.min(), ct_arr.max())


def populate_histogram(input_dir, mag, histogram):
    list_of_files = os.listdir(input_dir)

    dtype = np.dtype([('mag', float), ('ct', float)])

    for file_name in list_of_files:
        if file_name.startswith(mag):
            data = np.genfromtxt(os.path.join(input_dir, file_name), dtype=dtype)
            if 'mag' in histogram:
                np.testing.assert_array_equal(histogram['mag'], data['mag'])
                histogram['ct'] += data['ct']
            else:
                histogram['mag'] = data['mag']
                histogram['ct'] = data['ct']

    return None


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, default=None)
    parser.add_argument("--output_dir", type=str, default=None)

    args = parser.parse_args()
    if args.input_dir is None:
        raise RuntimeError("Must specify input_dir")
    if args.output_dir is None:
        raise RuntimeError("Must specify output_dir")
    if os.path.exists(args.output_dir) and not os.path.isdir(args.output_dir):
        raise RuntimeError("%s is not a dir" % args.output_dir)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    mag_list = ['u', 'g', 'r', 'i', 'z', 'y']

    color_color_fit_dir = os.path.join(args.input_dir, 'color_color_fit')
    color_color_input_dir = os.path.join(args.input_dir, 'color_color_input')

    for im in range(len(mag_list)-2):
        mag1 = mag_list[im]
        mag2 = mag_list[im+1]
        mag3 = mag_list[im+2]

        fig_name = os.path.join(args.output_dir, '%s_color_color_plots.png' % mag1)

        plt.figsize = (30, 30)

        xlabel = '%s-%s' % (mag1, mag2)
        ylabel = '%s-%s' % (mag2, mag3)

        grid = {}
        populate_grid(color_color_input_dir, mag1, grid)
        title = 'input'
        xbounds=(0,10)
        ybounds=(0,10)
        _xbounds, _ybounds, cbounds = plot_grid(grid, xlabel, ylabel, title, 1,
                                                xbounds=xbounds, ybounds=ybounds)

        grid = {}
        populate_grid(color_color_fit_dir, mag1, grid)
        title = 'fit'
        plot_grid(grid, xlabel, ylabel, title, 2,
                  xbounds=xbounds, ybounds=ybounds, cbounds=cbounds)

        plt.tight_layout()
        plt.savefig(fig_name)
        plt.close()


    fig_name = os.path.join(args.output_dir, 'magnitude_histograms.png')
    plt.figsize=(30,30)
    histogram_fit_dir = os.path.join(args.input_dir, 'fit_histogram')
    histogram_input_dir = os.path.join(args.input_dir, 'input_histogram')
    for i_fig, mag in enumerate(mag_list):
        fit_hist = {}
        populate_histogram(histogram_fit_dir, mag, fit_hist)
        input_hist = {}
        populate_histogram(histogram_input_dir, mag, input_hist)

        plt.subplot(3, 2, i_fig+1)
        area = (input_hist['ct']*(input_hist['mag'][1]-input_hist['mag'][0])).sum()
        hinput, = plt.plot(input_hist['mag'], input_hist['ct'], color='k')
        area = (fit_hist['ct']*(fit_hist['mag'][1]-fit_hist['mag'][0])).sum()
        hfit, = plt.plot(fit_hist['mag'], fit_hist['ct'], color='r')

        plt.xlabel('%s' % mag)
        plt.ylabel('N_stars')

        if i_fig == 0:
            plt.legend([hinput, hfit],
                       ['input', 'fit'], loc=0)

        print '%s fit: %d input: %d' % (mag, int(fit_hist['ct'].sum()), int(input_hist['ct'].sum()))

    plt.tight_layout()
    plt.savefig(fig_name)
    plt.close()

    fig_name = os.path.join(args.output_dir, 'input_vs_fit_mags.png')
    plt.figsize=(30,30)
    input_vs_fit_dir = os.path.join(args.input_dir, "input_vs_fit")
    for i_fig, mag in enumerate(mag_list):
        grid = {}
        populate_grid(input_vs_fit_dir, mag, grid)
        xbounds, ybounds, cbounds = plot_grid(grid, '%s_in' % mag, '%s_fit' % mag, None, i_fig+1,
                                              rows=3, cols=2)

        plt.plot(ybounds, ybounds, linewidth=0.5, color='r')

    plt.tight_layout()
    plt.savefig(fig_name)
    plt.close()
