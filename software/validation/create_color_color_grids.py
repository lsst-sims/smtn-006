from __future__ import with_statement
import argparse
import numpy as np
import os
import time

def _populate_grid(xx, yy, dd, grid):
    data_x = np.round(xx/dd).astype(int)
    data_y = np.round(yy/dd).astype(int)

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


def populate_mag_grid(data_in, mag, dmag, grid):
    fit_tag = '%snoatm' % mag
    input_tag = '%s_in' % mag
    if len(data_in)==0:
        return None
    _populate_grid(data_in[input_tag], data_in[fit_tag],
                   dmag, grid)

    return None


def populate_color_grid(data_in, mag1, mag2, mag3, dmag, grid):
    good_dexes = np.where(np.logical_and(data_in[mag1]>-98.0,
                          np.logical_and(data_in[mag2]>-98.0,
                                         data_in[mag3]>-98.0)))

    data = data_in[good_dexes]

    if len(data)==0:
        return None

    color1 = data[mag1]-data[mag2]
    color2 = data[mag2]-data[mag3]

    _populate_grid(color1, color2, dmag, grid)

    return None

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, default=None)
    parser.add_argument("--output_dir", type=str, default="plot_grids")
    parser.add_argument("--prefix", type=str, default=None)
    parser.add_argument("--min_colors", type=int, default=-1)
    parser.add_argument("--to_do_list", type=str, default=None)

    args = parser.parse_args()
    if args.input_dir is None:
        raise RuntimeError("Must specify input directory")
    if args.prefix is None:
        raise RuntimeError("Must specify prefix for output file")
    if args.to_do_list is None:
        raise RuntimeError("Must pass in a to_do_list")

    if not os.path.exists(args.output_dir) or not os.path.isdir(args.output_dir):
        raise RuntimeError("%s either does not exist or is not a dir" % args.output_dir)

    list_of_files = []
    with open(args.to_do_list, 'r') as input_file:
        lines = input_file.readlines()
        for line in lines:
            list_of_files.append(line.strip())

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

    fit_color_grids = {}
    input_color_grids = {}
    fit_input_mag_grids = {}
    mags = ['u', 'g', 'r', 'i', 'z', 'y']
    for ix in range(len(mags)-2):
        mag1 = mags[ix]
        fit_color_grids[mag1] = {}
        input_color_grids[mag1] = {}

    for mag1 in mags:
        fit_input_mag_grids[mag1] = {}

    mag_min = -10.0
    mag_max = 50.0
    fit_histograms = {}
    input_histograms = {}
    n_hist = int(np.round((mag_max-mag_min)/dmag))
    for mm in mags:
        fit_histograms[mm] = np.zeros(n_hist)
        input_histograms[mm] = np.zeros(n_hist)

    ct = 0
    ct_good = 0
    t_start = time.time()
    for file_name in list_of_files:
        if 'ebv_grid' not in file_name:
            continue

        full_name = os.path.join(args.input_dir, file_name)
        data = np.genfromtxt(full_name, dtype=dtype)

        good_dexes = np.where(np.logical_and(data['goodbit']==1, data['n_colors']>=args.min_colors))
        good_data = data[good_dexes]

        ct+=len(data)
        ct_good+=len(good_dexes[0])

        for mm in mags:
            populate_mag_grid(good_data, mm, dmag, fit_input_mag_grids[mm])

        for ix in range(len(mags)-2):
            mag1 = '%snoatm' % mags[ix]
            mag2 = '%snoatm' % mags[ix+1]
            mag3 = '%snoatm' % mags[ix+2]
            populate_color_grid(good_data, mag1, mag2, mag3, dmag, fit_color_grids[mags[ix]])

            mag1 = '%s_in' % mags[ix]
            mag2 = '%s_in' % mags[ix+1]
            mag3 = '%s_in' % mags[ix+2]
            populate_color_grid(good_data, mag1, mag2, mag3, dmag, input_color_grids[mags[ix]])

        for mm in mags:
            good_dexes = np.where(good_data['%snoatm' % mm]>-98.0)
            good_mags = good_data['%snoatm' % mm][good_dexes]
            ii_arr = (np.round((good_mags-mag_min)/dmag)).astype(int)
            unq, counts = np.unique(ii_arr, return_counts=True)
            for ii, cc in zip(unq, counts):
                fit_histograms[mm][ii] += cc

            good_dexes = np.where(good_data['%s_in' % mm]>-98.0)
            good_mags = good_data['%s_in' % mm][good_dexes]
            ii_arr = (np.round((good_mags-mag_min)/dmag)).astype(int)
            ii_arr = np.where(ii_arr>=0, ii_arr, 0)
            ii_arr = np.where(ii_arr<n_hist, ii_arr, n_hist-1)
            unq, counts = np.unique(ii_arr, return_counts=True)
            for ii, cc in zip(unq, counts):
                input_histograms[mm][ii] += cc

        print '%.3e %.3e time %.3e' % (ct,ct_good,time.time()-t_start)

    for mag1 in fit_color_grids:
        with open(os.path.join(args.output_dir, '%s_color_color_fit_%s.txt' % (args.prefix, mag1)), 'w') as out_file:
            grid = fit_color_grids[mag1]
            for xx in grid:
                for yy in grid[xx]:
                    out_file.write('%.3f %.3f %d\n' % (xx*dmag, yy*dmag, grid[xx][yy]))

    for mag1 in input_color_grids:
        with open(os.path.join(args.output_dir, '%s_color_color_input_%s.txt' % (args.prefix, mag1)), 'w') as out_file:
            grid = input_color_grids[mag1]
            for xx in grid:
                for yy in grid[xx]:
                    out_file.write('%.3f %.3f %d\n' % (xx*dmag, yy*dmag, grid[xx][yy]))

    for mm in mags:
        with open(os.path.join(args.output_dir, '%s_%s_fit_histogram.txt' % (args.prefix, mm)), 'w') as out_file:
            for ix in range(n_hist):
                out_file.write('%.3f %e\n' % (mag_min + ix*dmag, fit_histograms[mm][ix]))

        with open(os.path.join(args.output_dir, '%s_%s_input_histogram.txt' % (args.prefix, mm)), 'w') as out_file:
            for ix in range(n_hist):
                out_file.write('%.3f %e\n' % (mag_min + ix*dmag, input_histograms[mm][ix]))

        with open(os.path.join(args.output_dir, '%s_%s_input_vs_fit.txt' % (args.prefix, mm)), 'w') as out_file:
            grid = fit_input_mag_grids[mm]
            out_file.write('# input_mag fit_mag ct\n')
            for xx in grid:
                for yy in grid[xx]:
                    out_file.write('%.3f %.3f %d\n' % (xx*dmag, yy*dmag, grid[xx][yy]))
