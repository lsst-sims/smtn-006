"""
This script queries the CatSim database for all bright stars in
0 < ra < 20 and then returns density plots of declination versus
magnitude.
"""

import numpy as np
import matplotlib.pyplot as plt

from lsst.sims.catalogs.generation.db import DBObject

import time

if __name__ == "__main__":

    t_start = time.time()

    output_file_name = "mag_vs_dec.eps"

    db = DBObject(database='LSSTCATSIM', host='localhost',
                  port=51433, driver='mssql+pymssql')

    cmd = "SELECT catalogid, ra, decl, umag_noatm, gmag_noatm, rmag_noatm, "
    cmd += "imag_noatm, zmag_noatm, ymag_noatm, residual "
    cmd += "FROM bright_stars WHERE ra < 20.0 and ra > 0.0"

    query_dtype = np.dtype([('id', long),
                            ('ra', float), ('dec', float),
                            ('u', float), ('g', float), ('r', float),
                            ('i', float), ('z', float), ('y', float),
                            ('residual', float)])

    u_grid = {}
    g_grid = {}
    r_grid = {}
    i_grid = {}
    z_grid = {}
    y_grid = {}

    results = db.get_arbitrary_chunk_iterator(cmd, dtype=query_dtype,
                                              chunk_size=100000)

    dmag = 0.1
    ddec = 0.1

    total_stars = 0

    for chunk in results:
        print 'chunk ',time.time()-t_start,total_stars
        for line in chunk:
            total_stars += 1
            dec_int = int(round(line[2]/ddec))
            for grid_dict, mag_dex in \
            zip((u_grid, g_grid, r_grid, i_grid, z_grid, y_grid),
                (3, 4, 5, 6, 7, 8)):

                mag_int = int(round(line[mag_dex]/dmag))
                dict_key = '%s_%s' % (dec_int, mag_int)
                if dict_key in grid_dict:
                    grid_dict[dict_key][2] += 1
                else:
                    grid_dict[dict_key] = [ddec*dec_int, dmag*mag_int, 1]

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
            plt.xlabel('Dec')
            plt.ylabel('magnitude')

        if len(x_list)>2:
            xmax = x_list.max()
            xmin = x_list.min()
            ymax = y_list.max()
            ymin = y_list.min()

            plt.text(xmax-0.2*(xmax-xmin), ymax-0.1*(ymax-ymin), name, fontsize=15)


            dx = 0.1*(xmax-xmin)
            dy = 0.1*(ymax-ymin)

            xticks = np.arange(-90.0, 90.0, 5.0)
            xformat = ['%d' % xticks[ii] if ii%10==0 else '' for ii in range(len(xticks))]

            yticks = np.arange(np.floor(ymin),np.ceil(ymax),(np.ceil(ymax)-np.floor(ymin))*0.1)
            yformat = ['%.1f' % yticks[ii] if ii%3==0 else '' for ii in range(len(yticks))]

            plt.xticks(xticks, xformat, fontsize=15)
            plt.yticks(yticks, yformat, fontsize=15)

    plt.tight_layout()
    plt.savefig(output_file_name)
    plt.close()

    print 'that took ',time.time()-t_start
    print 'total stars ',total_stars
