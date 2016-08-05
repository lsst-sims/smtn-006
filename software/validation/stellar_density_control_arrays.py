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
import sys

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print "Need to pass in a to do list and a job dex"
        # job dex is an integer appended to the end of the output
        # files produced by this script so they don't collide with
        # the output files produced by another script.  This is done
        # so that you can run many jobs in parallel and get through
        # the whole catalog faster.


    to_do_list_file = sys.argv[1]
    job_dex = int(sys.argv[2])
    print to_do_list_file, job_dex

    t_start = time.time()

    # user-specified resolution.  Just make sure that it matches with
    # stellar_density_get_arrays.py and stellar_density_comparisons.py
    NSIDE = 64

    n_pix = hp.nside2npix(NSIDE)

    with open(to_do_list_file, "r") as input_file:
        input_file_list = input_file.readlines()

    mag_bounds = np.arange(9.0, 17.5, 0.5)

    grid_dict = {}
    for name in ('u', 'g', 'r', 'i', 'z', 'y'):
        grid_dict[name] = {}
        for mm in mag_bounds:
            grid_dict[name][mm] = np.zeros(n_pix)



    ra_max = None
    ra_min = None
    dec_max = None
    dec_min = None

    data_dtype = np.dtype([('id', long),
                           ('ra', float), ('dec', float),
                           ('mura', float), ('mudec', float),
                           ('b', float), ('v', float),
                           ('u', float), ('g', float),
                           ('r', float), ('i', float),
                           ('z', float), ('y', float),
                           ('j', float), ('h', float), ('k', float),
                           ('w1', float), ('w2', float),
                           ('w3', float), ('w4', float),
                           ('sst', float), ('flag', long)])

    total_stars = 0

    for file_name in input_file_list:
        file_name = file_name.strip()
        print file_name,total_stars,time.time()-t_start
        data = np.genfromtxt(file_name,
                             dtype=data_dtype, delimiter=',')

        total_stars += len(data['ra'])
        for name in ('u', 'g', 'r', 'i', 'z', 'y'):
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

                i_pix = np.unique(i_pix_raw)
                i_pix_add = np.array([len(np.where(i_pix_raw == ii)[0]) for ii in i_pix])
                grid_dict[name][bound][i_pix] += i_pix_add


    for name in ('u', 'g', 'r', 'i', 'z', 'y'):
      for ibound, mm in enumerate(mag_bounds):

            if ibound>0 and ibound<len(mag_bounds)-1:
                fig_title='%.2f <= %s < %.2f ' % (mag_bounds[ibound-1], name, mag_bounds[ibound])
            elif ibound==0:
                fig_title = '%s < %.2f' % (name, mag_bounds[ibound])
            else:
                fig_title = '%s > %.2f' % (name, mag_bounds[ibound])

            with open("healpix_control_grid_%s_%d_job%d.txt" % (name, ibound, job_dex), "w") as output_file:
                output_file.write("# %s\n" % fig_title)
                output_file.write("# nside %d\n" % NSIDE)

                for nn in grid_dict[name][mm]:
                    output_file.write("%d\n" % nn)

