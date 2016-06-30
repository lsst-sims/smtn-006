"""
This script crawls through sims_sed_library and computes the unnormalized,
un-extincted SDSS magnitudes and magNorms (magnitudes in the Imsim Bandpass)
of every stellar spectrum (main sequence, white dwarf, and M/L/T dwarf).

The results are output to the text file raw_sdss_mags.txt, which is meant
to be used by fit_bright_stars.cpp
"""

from __future__ import with_statement
import os
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Bandpass, Sed
import time


if __name__ == "__main__":

    # Load the SDSS bandpasses
    bp_dir = os.path.join(getPackageDir('throughputs'), 'sdss')
    bp_name_list = ['doi_u.dat', 'doi_g.dat', 'doi_r.dat', 'doi_i.dat',
                    'doi_z.dat']

    bp_list = []
    for name in bp_name_list:
        bp = Bandpass()
        bp.readThroughput(os.path.join(bp_dir, name))
        bp_list.append(bp)

    # load the imsim bandpass
    imsimBand = Bandpass()
    imsimBand.imsimBandpass()

    star_dir = os.path.join(getPackageDir('sims_sed_library'), 'starSED')
    sub_dir_list = ['kurucz', 'mlt', 'wDs']
    
    t_start = time.time()

    # loop through the kurucz, mlt, and wDs sub-directories of
    # sims_sed_library/starSED, calculating the required magnitudes
    # for each of the SED files found.
    with open('raw_sdss_mags.txt','w') as output_file:
        for sub_dir in sub_dir_list:
            sed_dir = os.path.join(star_dir, sub_dir)
            list_of_files = os.listdir(sed_dir)
            for file_name in list_of_files:
                output_file.write('%s ' % file_name)
                full_name = os.path.join(sed_dir, file_name)
                ss = Sed()
                ss.readSED_flambda(full_name)
                mm=ss.calcMag(imsimBand)
                output_file.write('%f ' % mm)

                for bp in bp_list:
                    mm = ss.calcMag(bp)
                    output_file.write('%f ' % mm)
                output_file.write('\n')
            print 'done with ',sub_dir,time.time()-t_start
