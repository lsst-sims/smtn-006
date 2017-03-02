from __future__ import with_statement
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import argparse
import numpy as np
import os
import time

def populate_grid(data_in, mag1, mag2, mag3, dmag, grid):
    good_dexes = np.where(np.logical_and(data_in[mag1]>-98.0,
                          np.logical_and(data_in[mag2]>-98.0,
                                         data_in[mag3]>-98.0)))

    data = data_in[good_dexes]

    if len(data)==0:
        return None

    color1 = data[mag1]-data[mag2]
    color2 = data[mag2]-data[mag3]

    data_x = np.round(color1/dmag).astype(int)
    data_y = np.round(color2/dmag).astype(int)

    xmax = data_x.max()
    ymax = data_y.max()
    xmin = data_x.min()
    ymin = data_y.min()
    factor = int(np.power(10.0, np.round(np.log10(ymax-ymin)+1.0)))
    dex_arr = (data_x-xmin)*factor + data_y-ymin

    unq, counts = np.unique(dex_arr, return_counts=True)

    x_unq = xmin + unq//factor
    y_unq = ymin + unq % factor

    for xx, yy, cc in zip(x_unq, y_unq, counts):
        if xx not in grid:
            grid[xx]= {}

        if yy not in grid[xx]:
            grid[xx][yy] = cc
        else:
            grid[xx][yy] += cc

    return None

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, default=None)
    parser.add_argument("--output_dir", type=str, default="plots")
    parser.add_argument("--prefix", type=str, default=None)

    args = parser.parse_args()
    if args.input_dir is None:
        raise RuntimeError("Must specify input directory")
    if args.prefix is None:
        raise RuntimeError("Must specify prefix for output file")

    dmag = 0.05

    dtype = np.dtype([('id', long), ('ra', float), ('dec', float),
                      ('pmra', float), ('pmdec', float),
                      ('lon', float), ('lat', float),
                      ('sedname', str, 300), ('magnorm', float),
                      ('flux_factor', float), ('ebv', float),
                      ('unoatm', float), ('gnoatm', float), ('rnoatm', float),
                      ('inoatm', float), ('znoatm', float), ('ynoatm', float),
                      ('u', float), ('g', float), ('r', float),
                      ('i', float), ('z', float), ('y', float),
                      ('b_in', float), ('v_in', float),
                      ('u_in', float), ('g_in', float), ('r_in', float),
                      ('i_in', float), ('z_in', float), ('y_in', float),
                      ('j_in', float), ('h_in', float), ('k_in', float),
                      ('w1_in', float), ('w2_in', float), ('w3_in', float),
                      ('w4_in', float), ('sst_in', float),
                      ('b4bit', int), ('2massbit', int), ('sdssbit', int),
                      ('ppmxlbit', int), ('ucac4bit', int), ('tychobit', int),
                      ('hipparcosbit', int), ('allwisebit', int),
                      ('levinesharabit', int), ('apassbit', int),
                      ('ps1bit', int), ('2massexbit', int),
                      ('astro0bit', int), ('astro1bit', int),
                      ('astro2bit', int), ('astro3bit', int),
                      ('photo0bit', int), ('photo1bit', int),
                      ('photo2bit', int), ('photo3bit', int),
                      ('clipbit', int), ('ppmxlresid_bit', int),
                      ('worrybit', int), ('gribit', int),
                      ('dirtybit', int), ('spikebit', int),
                      ('confusedbit', int), ('bvbit', int),
                      ('goodbit', int), ('norbit', int),
                      ('galaxbit', int), ('gvsbit', int),
                      ('residual', float), ('source', str, 300),
                      ('n_colors', int)])

    fit_grids = {}
    input_grids = {}
    mags = ['u', 'g', 'r', 'i', 'z', 'y']
    for ix in range(len(mags)-2):
        mag1 = mags[ix]
        fit_grids[mag1] = {}
        input_grids[mag1] = {}

    ct = 0
    ct_good = 0
    t_start = time.time()
    list_of_files = os.listdir(args.input_dir)
    for file_name in list_of_files:
        if 'ebv_grid' not in file_name:
            continue

        full_name = os.path.join(args.input_dir, file_name)
        data = np.genfromtxt(full_name, dtype=dtype)
        good_dexes = np.where(data['goodbit']==1)
        good_data = data[good_dexes]

        for ix in range(len(mags)-2):
            mag1 = '%snoatm' % mags[ix]
            mag2 = '%snoatm' % mags[ix+1]
            mag3 = '%snoatm' % mags[ix+2]
            populate_grid(good_data, mag1, mag2, mag3, dmag, fit_grids[mags[ix]])

            mag1 = '%s_in' % mags[ix]
            mag2 = '%s_in' % mags[ix+1]
            mag3 = '%s_in' % mags[ix+2]
            populate_grid(good_data, mag1, mag2, mag3, dmag, input_grids[mags[ix]])

        ct+=len(data)
        ct_good+=len(good_dexes[0])
        print '%.3e %.3e time %.3e' % (ct,ct_good,time.time()-t_start)

    for mag1 in fit_grids:
        with open(os.path.join(args.output_dir, '%s_%s_fit_grid_data.txt' % (args.prefix, mag1)), 'w') as out_file:
            grid = fit_grids[mag1]
            for xx in grid:
                for yy in grid[xx]:
                    out_file.write('%.6e %.6e %d\n' % (xx*dmag, yy*dmag, grid[xx][yy]))

    for mag1 in input_grids:
        with open(os.path.join(args.output_dir, '%s_%s_input_grid_data.txt' % (args.prefix, mag1)), 'w') as out_file:
            grid = input_grids[mag1]
            for xx in grid:
                for yy in grid[xx]:
                    out_file.write('%.6e %.6e %d\n' % (xx*dmag, yy*dmag, grid[xx][yy]))



    fig_name_list = [os.path.join(args.output_dir, '%s_fit_grid.png' % args.prefix),
                     os.path.join(args.output_dir, '%s_input_grid.png' % args.prefix)]
    grid_list = [fit_grids, input_grids]

    for fig_name, grid_set in zip(fig_name_list, grid_list):
        plt.figsize = (30, 30)
        for i_fig in range(len(mags)-2):
            plt.subplot(2, 2, i_fig+1)
            mag1 = mags[i_fig]
            mag2 = mags[i_fig+1]
            mag3 = mags[i_fig+2]

            grid =grid_set[mag1]

            x_arr = []
            y_arr = []
            ct_arr = []
            for xx in grid:
                for yy in grid[xx]:
                    x_arr.append(xx*dmag)
                    y_arr.append(yy*dmag)
                    ct_arr.append(grid[xx][yy])
            x_arr = np.array(x_arr)
            y_arr = np.array(y_arr)
            ct_arr = np.array(ct_arr)

            if len(ct_arr)==0:
                continue

            sorted_dexes = np.argsort(ct_arr)
            plt.scatter(x_arr[sorted_dexes], y_arr[sorted_dexes],
                        c=ct_arr[sorted_dexes],
                        edgecolor='',s=5,
                        cmap=plt.cm.gist_ncar,
                        norm=Normalize(vmin=ct_arr.min(), vmax=ct_arr.max()))

            cb = plt.colorbar()

            plt.xlabel('%s-%s' % (mag1, mag2))
            plt.ylabel('%s-%s' % (mag2, mag3))

        plt.tight_layout()
        plt.savefig(fig_name)
        plt.close()
