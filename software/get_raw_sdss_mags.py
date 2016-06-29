from __future__ import with_statement
import os
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Bandpass, Sed
import time


if __name__ == "__main__":

    bp_dir = os.path.join(getPackageDir('throughputs'), 'sdss')
    bp_name_list = ['doi_u.dat', 'doi_g.dat', 'doi_r.dat', 'doi_i.dat',
                    'doi_z.dat']

    bp_list = []
    for name in bp_name_list:
        bp = Bandpass()
        bp.readThroughput(os.path.join(bp_dir, name))
        bp_list.append(bp)

    imsimBand = Bandpass()
    imsimBand.imsimBandpass()

    star_dir = os.path.join(getPackageDir('sims_sed_library'), 'starSED')
    sub_dir_list = ['kurucz', 'mlt', 'wDs']
    
    t_start = time.time()
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
