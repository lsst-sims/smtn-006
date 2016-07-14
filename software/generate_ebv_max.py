import sys
import gzip
import numpy as np
from lsst.sims.photUtils import EBVbase

if __name__ == "__main__":

    input_name = sys.argv[1]

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
    ee = EBVbase()
    ebv = ee.calculateEbv(equatorialCoordinates=np.array([np.radians(data['ra']), np.radians(data['dec'])]))
    out_name = input_name.split('/')[-1]+'_ebv_max.txt'
    np.savetxt(out_name, ebv)
