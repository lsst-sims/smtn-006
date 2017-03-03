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

def plot_grid(grid, xlabel, ylabel, title, i_fig):
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

    plt.subplot(1, 2, i_fig)
    plt.scatter(x_arr, y_arr, c=cum_arr.astype(float)/ct_arr.sum(),
                edgecolor='',s=5,
                cmap=plt.cm.gist_ncar)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if i_fig==1:
        cb = plt.colorbar()


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
        plot_grid(grid, xlabel, ylabel, title, 1)

        grid = {}
        populate_grid(color_color_fit_dir, mag1, grid)
        title = 'fit'
        plot_grid(grid, xlabel, ylabel, title, 2)

        plt.tight_layout()
        plt.savefig(fig_name)
        plt.close()
