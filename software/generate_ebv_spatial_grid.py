from __future__ import with_statement
import numpy as np

from lsst.sims.photUtils import EBVbase
from lsst.sims.utils import sphericalFromCartesian

if __name__ == "__main__":

    rng = np.random.RandomState(99)
    n_grid = 100000
    ee = EBVbase()

    pt_list = rng.normal(0.0, scale=1.0, size=(n_grid, 3))
    lon_list, lat_list = sphericalFromCartesian(pt_list)
    ebv_list = ee.calculateEbv(equatorialCoordinates=np.array([lon_list, lat_list]))

    print ebv_list.min(),ebv_list.max()

    with open('ebv_spatial_grid.txt', 'w') as output_file:
        output_file.write('# RA Dec integrated E(B-V)\n')
        for lon, lat, ebv in zip(lon_list, lat_list, ebv_list):
            output_file.write('%.12f %.12f %.12f\n' % (np.degrees(lon), np.degrees(lat), ebv))

