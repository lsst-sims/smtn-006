"""
This script crawls through sims_sed_library and computes the unnormalized,
un-extincted SDSS magnitudes and magNorms (magnitudes in the Imsim Bandpass)
of every stellar spectrum (main sequence, white dwarf, and M/L/T dwarf).

The results are output to the text file ebv_independent_data.txt, which is
meant to be used by fit_bright_stars.cpp
"""

from __future__ import with_statement
import numpy as np
import os
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Bandpass, Sed
import time

def get_kurucz_phys(sed_name):
    """
    Read in the name of a kurucz SED file.  Return it's
    T_eff, metallicity, log(surface gravity)
    """
    if sed_name[1] == 'p':
        metallicity_sgn = 1.0
    elif sed_name[1] == 'm':
        metallicity_sgn = -1.0
    else:
        raise RuntimeError('Cannot parse metallicity sign of %s' % sed_name)

    new_name = sed_name.replace('.','_').split('_')

    metallicity = 0.1*metallicity_sgn*float(new_name[0][2:])

    teff = float(new_name[-2])

    logg = 0.1*np.float(new_name[3][1:])

    return teff, metallicity, logg


def get_wd_phys(sed_name):
    """
    Read in the name of a white dwarf SED,
    return its T_eff, metallicity (which we don't actually have),
    and log(surface gravity)
    """
    new_name = sed_name.replace('.','_').split('_')
    teff = float(new_name[-2])
    if new_name[1]!='He':
        logg = 0.1*float(new_name[2])
    else:
        logg = 0.1*float(new_name[3])

    return teff, -999.0, logg


def get_mlt_phys(sed_name):
    """
    Read in the name of an M/L/T dwarf SED and return
    its T_eff, metallicity, and log(surface gravity)
    """

    new_name = sed_name.replace('+','-').replace('a','-').split('-')

    logg_sgn_dex = len(new_name[0])

    if sed_name[logg_sgn_dex] == '-':
        logg_sgn = 1.0
    elif sed_name[logg_sgn_dex] == '+':
        logg_sgn = -1.0
    else:
        raise RuntimeError('Cannot get logg_sgn for %s' % sed_name)

    metallicity_sgn_dex = len(new_name[0]) + len(new_name[1]) + 1

    if sed_name[metallicity_sgn_dex] == '-':
        metallicity_sgn = -1.0
    elif sed_name[metallicity_sgn_dex] == '+':
        metallicity_sgn = 1.0
    else:
        raise RuntimeError('Cannot get metallicity_sgn for %s' % sed_name)

    teff = 100.0*float(new_name[0][3:])
    metallicity = metallicity_sgn*float(new_name[2])
    logg = logg_sgn*float(new_name[1])

    return teff, metallicity, logg


def get_physical_characteristics(sed_name, sub_dir):
    """
    Read in the name of an SED file.
    Return (in this order) Teff, metallicity (FeH), log(g)
    """
    if 'kurucz' in sub_dir:
        return get_kurucz_phys(sed_name)
    elif sub_dir == 'wDs':
        return get_wd_phys(sed_name)
    elif sub_dir == 'mlt':
        return get_mlt_phys(sed_name)
    else:
        raise RuntimeError('Do not know how to get '
                           'physical characteristics for '
                           'sub_dir %s' % sub_dir)



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

    #star_dir = os.path.join(getPackageDir('sims_sed_library'), 'starSED')
    sub_dir_list = ['kurucz', 'mlt', 'wDs']
    star_dir = os.path.join('/astro', 'users', 'danielsf',
                            'bright_stars', 'regen_kurucz')
    
    t_start = time.time()

    # loop through the kurucz, mlt, and wDs sub-directories of
    # sims_sed_library/starSED, calculating the required magnitudes
    # for each of the SED files found.
    with open('ebv_independent_data.txt','w') as output_file:
        for sub_dir in sub_dir_list:
            sed_dir = os.path.join(star_dir, sub_dir)
            list_of_files = os.listdir(sed_dir)
            for file_name in list_of_files:
                phys = get_physical_characteristics(file_name, sub_dir)
                output_file.write('%s %f %f %f ' % (file_name, phys[0], phys[1], phys[2]))
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
