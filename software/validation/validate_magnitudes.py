"""
This script reads in a sub-set of the raw csv files provided by Dave Monet,
queries the CatSim database for the stars that come from those files, and
then generates density plots of (input mag vs model mag - input mag) and
(model mag - input mag vs color residual)
"""


from __future__ import with_statement

import os
import numpy as np
from lsst.sims.catalogs.generation.db import DBObject

import matplotlib.pyplot as plt

import time

if __name__ == "__main__":

    output_dir = "magnitude_out_dir"

    monet_root_dir = os.path.join("/astro", "store", "pogo1", "monetBrightStars")
    monet_dir_list = []
    for ii in range(6):
        name = os.path.join(monet_root_dir, "Disc_%d" % ii)
        monet_dir_list.append(name)

    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

    data_dtype = np.dtype([('id', long),
                           ('ra', float),
                           ('dec', float),
                           ('mura', float),
                           ('mudec', float),
                           ('b', float), ('v', float),
                           ('u', float), ('g', float),
                           ('r', float), ('i', float),
                           ('z', float), ('y', float),
                           ('j', float), ('h', float),
                           ('k', float), ('w1', float),
                           ('w2', float), ('w3', float),
                           ('w4', float), ('sst', float),
                           ('flag', long)])

    query_dtype = np.dtype([('id', long),
                            ('u', float), ('g', float), ('r', float),
                            ('i', float), ('z', float), ('y', float),
                            ('residual', float)])

    dmag = 0.1
    u_grid = {}
    g_grid = {}
    r_grid = {}
    i_grid = {}
    z_grid = {}
    y_grid = {}

    u_resid_grid = {}
    g_resid_grid = {}
    r_resid_grid = {}
    i_resid_grid = {}
    z_resid_grid = {}
    y_resid_grid = {}

    total_stars = 0
    t_start = time.time()
    for monet_dir in monet_dir_list:
        list_of_files = os.listdir(monet_dir)
        for source_file in list_of_files:
            if source_file.endswith(".csv.gz"):
                monet_data_raw = np.genfromtxt(os.path.join(monet_dir, source_file), dtype=data_dtype,
                                               delimiter=',')

                total_stars += len(monet_data_raw['id'])

                print "read raw monet data"

                sorted_dexes = np.argsort(monet_data_raw['id'])

                monet_data = monet_data_raw[sorted_dexes]

                print "sorted monet data"

                cmd = "SELECT catalogid, umag_noatm, gmag_noatm, rmag_noatm, "
                cmd += "imag_noatm, zmag_noatm, ymag_noatm, residual "
                cmd += "FROM bright_stars WHERE source_file = '%s' " % source_file
                cmd += "ORDER BY catalogid"

                results = db.get_arbitrary_chunk_iterator(cmd, dtype=query_dtype, chunk_size=100000)

                ct = -1
                for chunk in results:
                    print "    chunk ",ct,time.time()-t_start
                    for line in chunk:
                        ct += 1
                        assert line[0] == monet_data['id'][ct]

                        for data_name, model_dex, grid_dict, resid_dict in \
                        zip(('u', 'g', 'r' ,'i', 'z', 'y'),
                            (1, 2, 3, 4, 5, 6), (u_grid, g_grid, r_grid, i_grid, z_grid, y_grid),
                            (u_resid_grid, g_resid_grid, r_resid_grid, i_resid_grid, z_resid_grid, y_resid_grid)):

                            if monet_data[data_name][ct]>-90.0:
                                data_int = int(round(monet_data[data_name][ct]/dmag))
                                residual_int = int(round((line[model_dex]-monet_data[data_name][ct])/dmag))
                                dict_tag = '%d_%d' % (data_int, residual_int)
                                if dict_tag in grid_dict:
                                    grid_dict[dict_tag][2] += 1
                                else:
                                    grid_dict[dict_tag] = [data_int*dmag, residual_int*dmag, 1]

                                color_int = int(round(line[7]/dmag))
                                dict_tag = '%d_%d' % (residual_int, color_int)
                                if dict_tag in resid_dict:
                                    resid_dict[dict_tag][2] += 1
                                else:
                                    resid_dict[dict_tag] = [residual_int*dmag, color_int*dmag, 1]


    for name, grid in zip(('u', 'g', 'r', 'i', 'z', 'y'),
                          (u_grid, g_grid, r_grid, i_grid, z_grid, y_grid)):

        with open(os.path.join(output_dir, "%s_mag_grid.txt" % name), "w") as output_file:

            output_file.write("# input mag, fit-input mag, ct\n")
            for k in grid:
                output_file.write("%e %e %d\n" % (grid[k][0], grid[k][1], grid[k][2]))


    for name, grid in zip(('u', 'g', 'r', 'i', 'z', 'y'),
                          (u_resid_grid, g_resid_grid, r_resid_grid, i_resid_grid, z_resid_grid, y_resid_grid)):

        with open(os.path.join(output_dir, "%s_fit_grid.txt" % name), "w") as output_file:

            output_file.write("# fit-input mag, color residual, ct\n")
            for k in grid:
                output_file.write("%e %e %d\n" % (grid[k][0], grid[k][1], grid[k][2]))


    print 'fitting took ',time.time()-t_start
    print "total stars ", total_stars
