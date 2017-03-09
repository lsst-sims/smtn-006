"""
This script crawls through sims_sed_library and computes the unnormalized
magnitudes of every stellar SED (main sequence, white dwarf, and M/L/T dwarf)
in every bandpass contained in the input bright star catalog at a grid
of E(B-V) values.  It produces a text file magnitude_grid.txt.  The columns
of that file are

sed file name
E(B-V)
2MASS J magnitude
2MASS H magnitude
2MASS Ks magnitude
WISE w1 magnitude
WISE w2 magnitude
WISE w3 magnitude
WISE w4 magnitude
Johnson U magnitude
Johnson B magnitude
Johnson V magnitude
Hipparcos Hp magnitude
Tycho B magnitude
Tycho V magnitude
PanStarrs g magnitude
PanStarrs r magnitude
PanStarrs i magnitude
PanStarrs z magnitude
PanStarrs y magnitude
SDSS u magnitude
SDSS g magnitude
SDSS r magnitude
SDSS i magnitude
SDSS z magnitude
LSST u magnitude
LSST g magnitude
LSST r magnitude
LSST i magnitude
LSST z magnitude
LSST y magnitude
"""

from __future__ import with_statement
import numpy as np
import os
import time

from lsst.utils import getPackageDir
from lsst.sims.photUtils import Sed, Bandpass

def get_bandpass_list(sub_dir, file_name_list):
    """
    Read in a sub directory of the throughputs package and
    a list of file names.  Return a list of Bandpasses
    read from sub_dir/file_name.
    """
    root_dir = getPackageDir('throughputs')
    bandpass_dir = os.path.join(root_dir, sub_dir)
    output_list = []
    dtype = np.dtype([('wv', np.float), ('sb', np.float)])
    for file_name in file_name_list:
        full_name = os.path.join(bandpass_dir, file_name)
        data = np.genfromtxt(full_name, dtype=dtype)
        bp = Bandpass(wavelen=data['wv'], sb=data['sb'],
                      wavelen_min=data['wv'].min(),
                      wavelen_max=data['wv'].max(),
                      wavelen_step=np.diff(data['wv']).min())
        output_list.append(bp)

    return output_list

def calc_magnitudes(ss, a_x, b_x, ebv, bandpass_list):
    """
    ss is an Sed
    a_x and b_x are the extinction law parameter grids
    ebv is E(B-V)
    bandpass_list is a list of Bandpasses

    return a list of magnitudes
    """
    wv = np.copy(ss.wavelen)
    flambda = np.copy(ss.flambda)

    ext_wav, ext_flambda = ss.addCCMDust(a_x, b_x, ebv=ebv,
                               wavelen=wv, flambda=flambda)

    ext_wav, ext_fnu = ss.flambdaTofnu(wavelen=ext_wav, flambda=ext_flambda)

    mag_list = [ss.calcMag(bp, wavelen=ext_wav, fnu=ext_fnu)
                for bp in bandpass_list]

    return mag_list



if __name__ == "__main__":

    # create grid of E(B-V) values
    ebv_list = np.arange(0.001, 0.01, 0.001)
    ebv_list = np.append(ebv_list,np.arange(0.01, 0.3, 0.01))
    ebv_list = np.append(ebv_list, np.arange(0.3, 7.0, 0.1))

    root_sed_dir = getPackageDir('sims_sed_library')
    sub_dir_list = ('kurucz', 'mlt', 'wDs')

    t_before_bp = time.time()

    # assemble a dict of lists of Bandpasses
    bp_name_dict = {}
    bp_list_dict = {}
    bp_name_dict['2MASS'] = ('2MASS_J.dat', '2MASS_H.dat', '2MASS_Ks.dat')
    bp_list_dict['2MASS'] = get_bandpass_list('2MASS',
                              bp_name_dict['2MASS'])

    bp_name_dict['WISE'] = ('WISE_w1.dat', 'WISE_w2.dat', 'WISE_w3.dat',
                              'WISE_w4.dat')
    bp_list_dict['WISE'] = get_bandpass_list('WISE', bp_name_dict['WISE'])

    bp_name_dict['LSST'] = ('hardware_u.dat', 'hardware_g.dat',
                              'hardware_r.dat', 'hardware_i.dat',
                              'hardware_z.dat', 'hardware_y.dat')

    bp_list_dict['LSST'] = get_bandpass_list('baseline', bp_name_dict['LSST'])

    bp_name_dict['LSSTatm'] = ('total_u.dat', 'total_g.dat', 'total_r.dat',
                               'total_i.dat', 'total_z.dat', 'total_y.dat')
    bp_list_dict['LSSTatm'] = get_bandpass_list('baseline', bp_name_dict['LSSTatm'])

    bp_name_dict['Hipparcos-Tycho'] = ('hipparcos/hipparcos_hp.dat',
                                        'tycho/tycho_B.dat',
                                        'tycho/tycho_V.dat')

    bp_list_dict['Hipparcos-Tycho'] = get_bandpass_list('',
                                       bp_name_dict['Hipparcos-Tycho'])

    bp_name_dict['Johnson'] = ('johnson_U.dat', 'johnson_B.dat',
                                  'johnson_V.dat')

    bp_list_dict['Johnson'] = get_bandpass_list('johnson',
                                 bp_name_dict['Johnson'])

    bp_name_dict['PanStarrs'] = ('panStarrs_g.dat', 'panStarrs_r.dat',
                                   'panStarrs_i.dat', 'panStarrs_z.dat',
                                   'panStarrs_y.dat')

    bp_list_dict['PanStarrs'] = get_bandpass_list('panStarrs',
                                   bp_name_dict['PanStarrs'])

    bp_name_dict['SDSS'] = ('doi_u.dat', 'doi_g.dat', 'doi_r.dat',
                                 'doi_i.dat', 'doi_z.dat')
    bp_list_dict['SDSS'] = get_bandpass_list('sdss', bp_name_dict['SDSS'])

    print 'getting bandpasses took %f' % (time.time()-t_before_bp)

    bp_tag_list = ('2MASS', 'WISE', 'Johnson', 'Hipparcos-Tycho',
        'PanStarrs', 'SDSS', 'LSST', 'LSSTatm')

    t_start = time.time()

    with open('magnitude_grid.txt', 'w') as output_file:

        # create header line for output file
        output_file.write('#sed_name E(B-V) ')
        for bp_tag in bp_tag_list:
            for bp_name in bp_name_dict[bp_tag]:
                output_file.write('%s ' % bp_name.split('/')[-1])
        output_file.write('\n')

        # loop through the sub directories of sims_sed_library
        for sub_dir in sub_dir_list:
            if not sub_dir.startswith('high_temp'):
                sed_dir = os.path.join(root_sed_dir, 'starSED', sub_dir)
            else:
                sed_dir = os.path.join('/astro', 'users', 'danielsf',
                                       'bright_stars', 'regen_kurucz',
                                       sub_dir)
            list_of_files = os.listdir(sed_dir)

            for file_name in list_of_files:

                full_name = os.path.join(sed_dir, file_name)
                ss = Sed()
                ss.readSED_flambda(full_name)

                # set up the parameter grid for the dust model on this SED
                a_x, b_x = ss.setupCCMab()

                # loop over the grid of E(B-V) values
                for ebv in ebv_list:
                    output_file.write('%s %f ' % (file_name.split('/')[-1], ebv))

                    # loop over the bandpasses, calculating magnitudes for this SED
                    # at this value of E(B-V)
                    for bp_tag in bp_tag_list:

                        bp_list = bp_list_dict[bp_tag]
                        mags = calc_magnitudes(ss, a_x, b_x, ebv, bp_list)
                        for mm in mags:
                            output_file.write('%f ' % mm)

                    output_file.write('\n')

    print 'getting the mags took %f' % (time.time()-t_start)
    print 'n_ebv %d' % len(ebv_list)
