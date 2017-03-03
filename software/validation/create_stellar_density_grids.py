"""
This file reads in a list of names of the raw csv catalogs provided
by Dave Monet, reads in those files, assembles them into HEALPIX maps
segregated by magnitude, and outputs those maps as text files to be
read in by another script and plotted.  These text files will literally
just be a list of number representing the number of stars in each HEALPIXel
(i.e. the first number will be the number of stars in the 1st HEALPIXel, the
second number will be the number of stars in the 2nd HEALPIXel, etc. with the
numbering determined by healpy's default scheme when calling ang2pix.  As of
this writing, healpy defaults to the RING indexing scheme).
"""

from __future__ import with_statement

import numpy as np
import healpy as hp
import os

import time

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--to_do", type=str, default=None)
    parser.add_argument("--suffix", type=str, default='')
    parser.add_argument("--output_dir", type=str, default=None)
    parser.add_argument("--input_dir", type=str, default=None)
    parser.add_argument("--dmag", type=float, default=1.0)
    parser.add_argument("--min_colors", type=int, default=-1)

    args = parser.parse_args()
    if args.to_do is None:
        raise RuntimeError("Need to pass in a to do list")
    if args.output_dir is None:
        raise RuntimeError("Need to specify output_dir")
    if args.input_dir is None:
        raise RuntimeError("Need to specify input_dir")

    if os.path.exists(args.output_dir) and not os.path.isdir(args.output_dir):
        raise RuntimeError("%s is not a dir" % args.output_dir)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    fit_maps_dir = os.path.join(args.output_dir, "fit_healpix_maps")
    if not os.path.exists(fit_maps_dir):
        os.mkdir(fit_maps_dir)

    input_maps_dir = os.path.join(args.output_dir, "input_healpix_maps")
    if not os.path.exists(input_maps_dir):
        os.mkdir(input_maps_dir)

    mag_bounds = np.arange(9.0, 17.5, args.dmag)

    to_do_list_file = args.to_do

    t_start = time.time()

    # user-specified resolution.  Just make sure that it matches with
    # stellar_density_get_arrays.py and stellar_density_comparisons.py
    NSIDE = 64

    n_pix = hp.nside2npix(NSIDE)

    with open(to_do_list_file, "r") as input_file:
        input_file_list = input_file.readlines()

    mag_name_list = ['u_in', 'g_in', 'r_in', 'i_in', 'z_in', 'y_in',
                     'unoatm', 'gnoatm', 'rnoatm', 'inoatm', 'znoatm', 'ynoatm']

    grid_dict = {}
    for name in mag_name_list:
        grid_dict[name] = {}
        for mm in mag_bounds:
            grid_dict[name][mm] = np.zeros(n_pix)

    ra_max = None
    ra_min = None
    dec_max = None
    dec_min = None

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

    total_stars = 0

    for file_name in input_file_list:
        file_name = file_name.strip()
        print file_name,total_stars,time.time()-t_start
        data = np.genfromtxt(os.path.join(args.input_dir, file_name),
                             dtype=dtype, delimiter=' ')

        good_dexes = np.where(np.logical_and(data['n_colors']>=args.min_colors,
                                             data['goodbit']==1))

        data = data[good_dexes]

        total_stars += len(data['ra'])
        for name in mag_name_list:
            for ix, bound in enumerate(mag_bounds):
                if ix > 0 and ix<len(mag_bounds)-1:
                    min_mag = mag_bounds[ix-1]
                    max_mag = mag_bounds[ix]
                elif ix==0:
                    min_mag = -90.0
                    max_mag = mag_bounds[ix]
                else:
                    min_mag = mag_bounds[ix]
                    max_mag = 1000000.0


                valid_stars = np.where(np.logical_and(data[name]<max_mag, data[name]>=min_mag))
                i_pix_raw = hp.ang2pix(NSIDE, np.radians(90.0-data['dec'][valid_stars]),
                                       np.radians(data['ra'][valid_stars]))

                i_pix, counts = np.unique(i_pix_raw, return_counts=True)
                grid_dict[name][bound][i_pix] += counts

    for sub_dir in (fit_maps_dir, input_maps_dir):
        with open(os.path.join(sub_dir, 'number_of_bounds.txt'), 'w') as out_file:
            out_file.write('%d\n' % len(mag_bounds))

    for name in mag_name_list:
        if 'noatm' in name:
            sub_dir = fit_maps_dir
            print name,sub_dir
        elif '_in' in name:
            sub_dir = input_maps_dir
            print name,sub_dir
        file_name = name.replace('_in','').replace('noatm','')

        for ibound, mm in enumerate(mag_bounds):
            if ibound>0 and ibound<len(mag_bounds)-1:
                fig_title='%.2f <= %s < %.2f ' % (mag_bounds[ibound-1], file_name, mag_bounds[ibound])
            elif ibound==0:
                fig_title = '%s < %.2f' % (file_name, mag_bounds[ibound])
            else:
                fig_title = '%s > %.2f' % (file_name, mag_bounds[ibound])

            with open(os.path.join(sub_dir, "%s_%d%s.txt" %
                      (file_name, ibound, args.suffix)), "w") as output_file:

                output_file.write("# %s\n" % fig_title)
                output_file.write("# nside %d\n" % NSIDE)

                for nn in grid_dict[name][mm]:
                    output_file.write("%d\n" % nn)

