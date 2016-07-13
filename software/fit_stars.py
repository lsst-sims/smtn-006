from __future__ import with_statement
import sys
import numpy as np
import gzip
from lsst.sims.photUtils import EBVbase
from lsst.sims.utils import galacticFromEquatorial

import time

_hexadec_places = 8

def load_mag_grid(file_name):

    dtype = np.dtype([('sed_name', str, 300),
                      ('ebv', np.float),
                      ('2MASS_J', np.float), ('2MASS_H', np.float),
                      ('2MASS_Ks', np.float),
                      ('WISE_w1', np.float), ('WISE_w2', np.float),
                      ('WISE_w3', np.float), ('WISE_w4', np.float),
                      ('johnson_U', np.float), ('johnson_B', np.float),
                      ('johnson_V', np.float),
                      ('hipparcos_hp', np.float), ('tycho_B', np.float),
                      ('tycho_V', np.float),
                      ('panStarrs_g', np.float), ('panStarrs_r', np.float),
                      ('panStarrs_i', np.float), ('panStarrs_z', np.float),
                      ('panStarrs_y', np.float),
                      ('sdss_u', np.float), ('sdss_g', np.float),
                      ('sdss_r', np.float), ('sdss_i', np.float),
                      ('sdss_z', np.float),
                      ('lsst_noatm_u', np.float), ('lsst_noatm_g', np.float),
                      ('lsst_noatm_r', np.float), ('lsst_noatm_i', np.float),
                      ('lsst_noatm_z', np.float), ('lsst_noatm_y', np.float),
                      ('lsst_u', np.float), ('lsst_g', np.float),
                      ('lsst_r', np.float), ('lsst_i', np.float),
                      ('lsst_z', np.float), ('lsst_y', np.float)])

    sed_mag_grid = np.genfromtxt(file_name, dtype=dtype)
    sorted_dexes = np.argsort(sed_mag_grid['ebv'])

    ebv_sorted = sed_mag_grid['ebv'][sorted_dexes]
    ebv_cut_vals = []
    ebv_cut_dexes = []
    for ix, val in enumerate(ebv_sorted):
        if ix>0:
            if val-ebv_sorted[ix-1]>1.0e-6:
                ebv_cut_vals.append(ebv_sorted[ix-1])
                ebv_cut_dexes.append(ix-1)

    ebv_cut_vals = np.array(ebv_cut_vals)
    ebv_cut_dexes = np.array(ebv_cut_dexes)

    mag_name_list = ('2MASS_J', '2MASS_H', '2MASS_Ks', 'WISE_w1', 'WISE_w2',
                     'WISE_w3', 'WISE_w4', 'johnson_U', 'johnson_B',
                     'johnson_V', 'hipparcos_hp', 'tycho_B',
                      'tycho_V', 'panStarrs_g', 'panStarrs_r', 'panStarrs_i',
                      'panStarrs_z', 'panStarrs_y', 'sdss_u','sdss_g','sdss_r',
                      'sdss_i', 'sdss_z', 'lsst_noatm_u', 'lsst_noatm_g', 'lsst_noatm_r',
                      'lsst_noatm_i', 'lsst_noatm_z','lsst_noatm_y', 'lsst_u', 'lsst_g',
                      'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y')

    return (sed_mag_grid['sed_name'][sorted_dexes], ebv_sorted,
            np.array([sed_mag_grid[nn][sorted_dexes] for nn in mag_name_list]),
            ebv_cut_vals, ebv_cut_dexes)

    return np.array(sed_mag_grid_raw[sorted_dexes], dtype=dtype)


def load_ebv_independent_data(file_name):
    dtype = np.dtype([('sed_name', str, 300),
                      ('teff', np.float), ('feh', np.float), ('logg', np.float),
                      ('magnorm', np.float),
                      ('sdss_u', np.float), ('sdss_g', np.float),
                      ('sdss_r', np.float), ('sdss_i', np.float),
                      ('sdss_z', np.float)
                     ])

    data = np.genfromtxt(file_name, dtype=dtype)
    dex_dict = {}
    for ix, name in enumerate(data['sed_name']):
        dex_dict[name] = ix
    return data, dex_dict

def convert_catalog_line(line):
    vv = line.split(',')

    output_dict = {}
    output_dict['star_id'] = long(vv[0])
    output_dict['ra'] = np.float(vv[1])
    output_dict['dec'] = np.float(vv[2])
    output_dict['mura'] = np.float(vv[3])
    output_dict['mudec'] = np.float(vv[4])
    output_dict['flag'] = long(vv[21])
    output_dict['mags'] = np.array([np.float(nn) for nn in vv[5:20]])

    return output_dict


def twos_complement(ii_in):

    binary_places = 4*_hexadec_places
    binary_raw = np.zeros(binary_places, np.int)
    local_term = long(2**binary_places)
    remainder = long(-1*ii_in)
    for dex in range(binary_places-1, -1, -1):
        local_term/=2
        val=remainder/local_term
        if(val>1):
            raise RuntimeError("Failure with binary %lld %lld %d" %
                               (remainder,local_term,val))
        binary_raw[dex]=val
        remainder-=val*local_term


    binary = [1 if nn==0 else 0 for nn in binary_raw]

    if binary[0]==0:
        binary[0]=1
    else:
        binary[0]=0
        for i in range(1, binary_places, 1):
            if binary[i]==1:
                binary[i]=0
            else:
                binary[i]=1
                break

    remainder=long(0)
    local_term=long(1)
    for dex in range(binary_places):
        if dex>0:
            local_term*=2
        remainder += binary[dex]*local_term
    return remainder


def convert_to_hexadecimal(ii_in):

    if ii_in<0:
        ii = twos_complement(ii_in)
    else:
        ii = ii_in

    output = np.zeros(_hexadec_places, dtype=int)
    remainder = long(ii)
    local_term = long(16**_hexadec_places)
    for dex in range(_hexadec_places-1, -1, -1):
        local_term/=16
        val=remainder/local_term
        if val>=16:
            raise RuntimeError("Cannot convertd %lld to hexadecimal" % ii_in)

        output[dex]=val
        remainder-=val*local_term

    return output


def assemble_grid(flag, star_mags, raw_grid):

    hdec = convert_to_hexadecimal(flag)
    src_flag = hdec[4]

    name_grid = ['2MASS_J','2MASS_Ks', '2MASS_H', 'johnson_B',
                'johnson_V', 'sdss_u', 'sdss_g', 'sdss_r',
                'sdss_i', 'sdss_z', 'panStarrs_y', 'WISE_w1',
                'WISE_w2', 'WISE_w3', 'WISE_w4']

    dex_grid = [8, 9, 18, 19, 20, 21, 22, 17, 0, 1, 2, 3, 4, 5, 6]

    if src_flag in (0, 1, 11):
        name_grid[3] = 13
        name_grid[4] = 14
        name_grid[5] = 15
        name_grid[6] = 16
    elif src_flag in (6, 7):
        name_grid[0] = 11
        name_grid[1] = 12

    for ix in range(len(dex_grid)-1,-1,-1):
        if star_mags[ix]<-90.0:
            dex_grid.pop(ix)

    return (raw_grid[dex_grid].transpose(),
            np.array([mm for mm in star_mags if mm > -90.0]))



def fit_star(star_mags, grid_in, ebv_max, ebv_cut_vals, ebv_cut_dexes):

    if len(star_mags) != grid_in.shape[1]:
        raise RuntimeError("%d starmags; %d mags in grid" % (len(star_mags), grid_in.shape[1]))

    if len(ebv_grid) != grid_in.shape[0]:
        raise RuntimeError("%d ebv values; %d SEDs" % (len(ebv_grid), grid_in.shape[0]))

    #valid_ebv_dexes = np.where(ebv_grid <= ebv_max)
    #valid_ebv = ebv_grid[valid_ebv_dexes]
    #valid_names = name_grid[valid_ebv_dexes]

    _valid_dexes = np.where(ebv_cut_vals<=ebv_max+1.0e-6)
    cut_dex = ebv_cut_dexes[_valid_dexes[0][-1]]+1

    grid = grid_in[:cut_dex].transpose()
    #print grid.shape
    #print len(valid_ebv_dexes[0]),valid_ebv_dexes[0][0],valid_ebv_dexes[0][-1]
    #print valid_ebv_dexes[0]

    color_grid = (grid[1:]-grid[:-1]).transpose()

    #star_colors = np.array([star_mags[ix+1] - star_mags[ix] for ix in range(len(star_mags)-1)])
    star_colors = star_mags[1:]-star_mags[:-1]

    dex = np.argmin(np.power(star_colors-color_grid, 2).sum(axis=1))
    rms = np.sqrt(np.power(star_colors-color_grid[dex],2).sum())/len(star_colors)

    offset = (star_mags - grid_in[dex]).sum()/len(star_mags)
    return dex, rms, offset



if __name__ == "__main__":


    print "testing hexadec"
    test = convert_to_hexadecimal(431)
    control = np.zeros(8, np.int)
    control[0]=15
    control[1]=10
    control[2]=1
    np.testing.assert_array_equal(test, control)

    test = convert_to_hexadecimal(218108420)
    control = np.zeros(8,np.int)
    control[0]=4
    control[2]=2
    control[3]=1
    control[6]=13
    np.testing.assert_array_equal(test, control)

    test = convert_to_hexadecimal(-53)
    control = np.ones(8,np.int)
    control*=15
    control[0]=11
    control[1]=12
    np.testing.assert_array_equal(test,control)

    print "done testing hexadec"

    if len(sys.argv)<4:
        print "Call signature:"
        print "python fit_stars.py sed_mag_grid ebv_independent to_do"
        exit()

    t_start = time.time()

    sed_mag_grid_name = sys.argv[1]
    ebv_independent_name = sys.argv[2]
    to_do_list = sys.argv[3]

    name_grid, ebv_grid, mag_grid, ebv_cut_vals, ebv_cut_dexes = load_mag_grid(sed_mag_grid_name)

    # code to test that magnitude grid, etc is sorted correctly
    """
    tp_grid = mag_grid.transpose()
    with open(sed_mag_grid_name, 'r') as control_file:
        control_lines = control_file.readlines()
        for dex in range(0, len(name_grid), 17777):
            test_line = "%s %f" % (name_grid[dex], ebv_grid[dex])
            for tt in tp_grid[dex]:
                test_line += " %f" % tt
            test_line += ' \n'
            try:
                assert test_line in control_lines
            except:
                print 'could not find'
                print test_line
                raise
    """

    ebv_independent_data, ebv_independent_dexes = load_ebv_independent_data(ebv_independent_name)

    ebv_calculator = EBVbase()

    print "overhead took ",time.time()-t_start
    print mag_grid.shape

    out_name = "test_dir/test_output.txt"

    t_start = time.time()
    star_ct = 0
    with open(to_do_list, "r") as input_file_list:
        for file_name in input_file_list:
            file_name=file_name.strip()
            record_name = file_name.split('/')[-1]
            with gzip.open(file_name, "r") as catalog_file:
                with open(out_name, "w") as output_file:
                    output_file.write("#dummy_htmid star_id ra dec mura mudec lon lat ");
                    output_file.write("sed magnorm E(B-V) Teff [Fe/H] log(g) ");
                    output_file.write("lsst_u_noatm lsst_g_noatm lsst_r_noatm lsst_i_noatm lsst_z_noatm lsst_y_noatm ");
                    output_file.write("lsst_u_atm lsst_g_atm lsst_r_atm lsst_i_atm lsst_z_atm lsst_y_atm ");
                    output_file.write("sdss_u(ext) sdss_g(ext) sdss_r(ext) sdss_i(ext) sdss_z(ext) ");
                    output_file.write("sdss_u(raw) sdss_g(raw) sdss_r(raw) sdss_i(raw) sdss_z(raw) ");
                    output_file.write("color_residual file_name\n");
                    for line in catalog_file:
                        star_ct += 1
                        star_dict = convert_catalog_line(line)
                        lon, lat = galacticFromEquatorial(star_dict['ra'], star_dict['dec'])
                        ebv_max = ebv_calculator.calculateEbv(equatorialCoordinates=np.array([[np.radians(star_dict['ra'])],
                                                                                              [np.radians(star_dict['dec'])]]))

                        fit_grid, star_mags = assemble_grid(star_dict['flag'], star_dict['mags'], mag_grid)

                        dex,  rms, offset = fit_star(star_mags, fit_grid,
                                                     ebv_max, ebv_cut_vals, ebv_cut_dexes)

                        name = name_grid[dex]
                        ebv = ebv_grid[dex]
                        ind_dex = ebv_independent_dexes[name]

                        output_file.write("0 %ld %le %le %le %le %le %le " %
                                          (star_dict['star_id'], star_dict['ra'], star_dict['dec'],
                                           star_dict['mura']*1000.0, star_dict['mudec']*1000.0,
                                           lon,lat))

                        output_file.write("%s %le %le %le %le %le " %
                                          (name, ebv_independent_data['magnorm'][ind_dex] + offset, ebv,
                                           ebv_independent_data['teff'][ind_dex],
                                           ebv_independent_data['feh'][ind_dex],
                                           ebv_independent_data['logg'][ind_dex]))

                        # lsst noatm
                        output_file.write("%le %le %le %le %le %le " %
                                          (mag_grid[23][dex]+offset, mag_grid[24][dex]+offset, mag_grid[25][dex]+offset,
                                           mag_grid[26][dex]+offset, mag_grid[27][dex]+offset, mag_grid[28][dex]+offset))

                        # lsst atm
                        output_file.write("%le %le %le %le %le %le " %
                                          (mag_grid[29][dex]+offset, mag_grid[30][dex]+offset, mag_grid[31][dex]+offset,
                                           mag_grid[32][dex]+offset, mag_grid[33][dex]+offset, mag_grid[34][dex]+offset))


                        # sdss ext
                        output_file.write("%le %le %le %le %le " %
                                          (mag_grid[18][dex]+offset, mag_grid[19][dex]+offset, mag_grid[20][dex]+offset,
                                           mag_grid[21][dex]+offset, mag_grid[22][dex]+offset))

                        # sdss noext
                        output_file.write("%le %le %le %le %le " %
                                          (ebv_independent_data['sdss_u'][ind_dex] + offset,
                                           ebv_independent_data['sdss_g'][ind_dex] + offset,
                                           ebv_independent_data['sdss_r'][ind_dex] + offset,
                                           ebv_independent_data['sdss_i'][ind_dex] + offset,
                                           ebv_independent_data['sdss_z'][ind_dex] + offset))

                        output_file.write("%le %s \n" % (rms, record_name))


    print 'actual work took ',time.time()-t_start,' to do ',star_ct
