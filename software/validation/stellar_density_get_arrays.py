"""
This script queries the CatSim database for bright stars, separates them into
HEALPIX maps based on magnitudes, and outputs those HEALPIX maps as text files
to be read in and plotted by another script.  These text files will literally
just be a list of number representing the number of stars in each HEALPIXel
(i.e. the first number will be the number of stars in the 1st HEALPIXel, the
second number will be the number of stars in the 2nd HEALPIXel, etc. with the
numbering determined by healpy's default scheme when calling ang2pix.  As of
this writing, healpy defaults to the RING indexing scheme).
"""

from __future__ import with_statement

import numpy as np
import healpy as hp
from lsst.sims.catalogs.generation.db import DBObject

import time

if __name__ == "__main__":

    t_start = time.time()

    # user-specified resolution.  Just make sure that it matches with
    # stellar_density_comparisons.py and stellar_density_control_arrays.py
    NSIDE = 64

    n_pix = hp.nside2npix(NSIDE)

    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

    mag_bounds = np.arange(9.0, 17.5, 0.5)

    grid_dict = {}
    for name in ('u', 'g', 'r', 'i', 'z', 'y'):
        grid_dict[name] = {}
        for mm in mag_bounds:
            grid_dict[name][mm] = np.zeros(n_pix)

    print 'n_pix ',n_pix


    ra_max = None
    ra_min = None
    dec_max = None
    dec_min = None

    cmd = "SELECT ra, decl, umag_noatm, gmag_noatm, rmag_noatm, imag_noatm, "
    cmd += "zmag_noatm, ymag_noatm "
    cmd += "FROM bright_stars"

    query_dtype = np.dtype([('ra', float), ('dec', float),
                            ('u', float), ('g', float), ('r', float),
                            ('i', float), ('z', float), ('y', float)])

    results = db.get_arbitrary_chunk_iterator(cmd, dtype=query_dtype, chunk_size=100000)

    total_stars = 0

    for chunk in results:
        print "chunk ",total_stars,time.time()-t_start
        total_stars += len(chunk['ra'])
        for name in ('u', 'g', 'r', 'i', 'z', 'y'):
            for ix, bound in enumerate(mag_bounds):
                if ix > 0 and ix<len(mag_bounds)-1:
                    min_mag = mag_bounds[ix-1]
                    max_mag = mag_bounds[ix]
                elif ix==0:
                    min_mag = -99.0
                    max_mag = mag_bounds[ix]
                else:
                    min_mag = mag_bounds[ix]
                    max_mag = 1000000.0


                valid_stars = np.where(np.logical_and(chunk[name]<max_mag, chunk[name]>=min_mag))
                i_pix_raw = hp.ang2pix(NSIDE, np.radians(90.0-chunk['dec'][valid_stars]),
                                       np.radians(chunk['ra'][valid_stars]))

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

            with open("healpix_grid_%s_%d.txt" % (name, ibound), "w") as output_file:
                output_file.write("# %s\n" % fig_title)
                output_file.write("# nside %d\n" % NSIDE)

                for nn in grid_dict[name][mm]:
                    output_file.write("%d\n" % nn)

