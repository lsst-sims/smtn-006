"""
This script is called by fit_bright_stars.cpp to calculate the integrated
E(B-V) dust extinction along a line of sight.  fit_bright_stars.cpp then uses
that value as a cutoff when searching the E(B-V) grid for best fit values for
particular stars.

The call signature of the script is

python generate_ebv_max.py input_file output_file

where input_file is one of the original .csv.gz files as provided by Dave
Monet and output_file is a text file where, for each star in input_file, the
maximum E(B-V) value will be listed.
"""

import sys
import gzip
import numpy as np
from lsst.sims.photUtils import EBVbase

if __name__ == "__main__":

    input_name = sys.argv[1]
    output_name = sys.argv[2]

    dtype = np.dtype([('id', long),
                      ('ra', np.float), ('dec', np.float),
                      ('mura', np.float), ('mudec', np.float),
                      ('b', np.float), ('v', np.float),
                      ('u', np.float), ('g', np.float), ('r', np.float),
                      ('i', np.float), ('z', np.float), ('y', np.float),
                      ('j', np.float), ('h', np.float), ('k', np.float),
                      ('w1', np.float), ('w2', np.float),
                      ('w3', np.float), ('w4', np.float),
                      ('sst', np.float), ('flag', long)])

    data = np.genfromtxt(input_name, dtype=dtype, delimiter=',')
    ee = EBVbase()   # an instantiation of EBVbase from which to do calculations
    ebv = ee.calculateEbv(equatorialCoordinates=np.array([np.radians(data['ra']), np.radians(data['dec'])]))
    np.savetxt(output_name, ebv)
