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

    mag_file_name = "magnitude_fit_plots.eps"
    color_file_name = "goodness_of_fit_plots.eps"

    monet_dir = "monetData"  # the directory where the csv files are stored
    list_of_monet = os.listdir(monet_dir)
    list_of_monet.sort()

    db = DBObject(database='LSSTCATSIM', host='localhost',
                  port=51433, driver='mssql+pymssql')

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

    for source_file in list_of_monet:
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

        t_start = time.time()
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

        print 'fitting took ',time.time()-t_start

    plt.figsize = (30, 30)

    for i_fig, (grid_dict, name) in \
    enumerate(zip((u_grid, g_grid, r_grid, i_grid, z_grid, y_grid),
                  ('u', 'g', 'r', 'i', 'z', 'y'))):

        plt.subplot(3, 2, i_fig+1)

        x_list = []
        y_list = []
        color_list = []
        for key in grid_dict:
            x_list.append(grid_dict[key][0])
            y_list.append(grid_dict[key][1])
            color_list.append(grid_dict[key][2])

        x_list = np.array(x_list)
        y_list = np.array(y_list)
        color_list = np.array(color_list)

        color_dex = np.argsort(color_list)

        plt.scatter(x_list[color_dex], y_list[color_dex],
                   c=np.log10(color_list[color_dex]), s=5,
                   cmap=plt.cm.gist_ncar, edgecolor='')

        plt.colorbar()

        if i_fig==0:
            plt.xlabel('$m_{input}$')
            plt.ylabel('$m_{fit}-m_{input}$')

        if len(x_list)>2:
            xmax = x_list.max()
            xmin = x_list.min()
            ymax = y_list.max()
            ymin = y_list.min()

            plt.text(xmax-0.2*(xmax-xmin), ymax-0.1*(ymax-ymin), name, fontsize=15)

    plt.tight_layout()
    plt.savefig(mag_file_name)
    plt.close()


###########now the color residuals versus mag residuals
    plt.figsize = (30, 30)

    for i_fig, (grid_dict, name) in \
    enumerate(zip((u_resid_grid, g_resid_grid, r_resid_grid,
                   i_resid_grid, z_resid_grid, y_resid_grid),
                  ('u', 'g', 'r', 'i', 'z', 'y'))):

        plt.subplot(3, 2, i_fig+1)

        x_list = []
        y_list = []
        color_list = []
        for key in grid_dict:
            x_list.append(grid_dict[key][0])
            y_list.append(grid_dict[key][1])
            color_list.append(grid_dict[key][2])

        x_list = np.array(x_list)
        y_list = np.array(y_list)
        color_list = np.array(color_list)

        color_dex = np.argsort(color_list)

        plt.scatter(x_list[color_dex], y_list[color_dex],
                   c=np.log10(color_list[color_dex]), s=5,
                   cmap=plt.cm.gist_ncar, edgecolor='')

        plt.colorbar()

        if len(x_list)>2:
            xmax = x_list.max()
            xmin = x_list.min()
            ymax = y_list.max()
            ymin = y_list.min()

            plt.text(xmax-0.2*(xmax-xmin), ymax-0.1*(ymax-ymin), name, fontsize=15)

        if i_fig == 0:
            plt.xlabel('$m_{fit}-m_{input}$')
            plt.ylabel('color residual')

    plt.tight_layout()
    plt.savefig(color_file_name)

    print "total stars ", total_stars
