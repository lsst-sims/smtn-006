import matplotlib
matplotlib.use('Agg')

import os
import numpy as np

import matplotlib.pyplot as plt

def get_cumulative(in_mag, in_ct):

    in_tot = in_ct.sum()

    mag_resid = np.unique(np.abs(in_mag))
    mag_ct = np.array([in_ct[np.where(np.abs(in_mag) == mm)].sum() for mm in mag_resid])
    sorted_dex = np.argsort(mag_resid)
    mag_resid = mag_resid[sorted_dex]
    mag_ct = mag_ct[sorted_dex]

    cum_ct = [mag_ct[ii:].sum() for ii in range(len(mag_resid))]

    tot = mag_ct.sum()
    assert tot == cum_ct[0]
    assert tot == in_tot
    return mag_resid, np.array(cum_ct)


if __name__ == "__main__":

    data_dir = "mag_valid_grids"
    output_dir = "mag_figs"

    color_cuts = (3.0, 1.0, 0.75, 0.5)
    color_colors = ('r', 'b', 'm', 'g')

    name_list = ('u', 'g', 'r', 'i', 'z', 'y')

    mag_dtype = np.dtype([('input', float), ('mag_resid', float), ('ct', int)])
    fit_dtype = np.dtype([('mag_resid', float),
                          ('color_resid', float),
                          ('ct', int)])


    mag_data = {}
    fit_data = {}
    mag_resid_dist = {}

    ###read in data
    for name in name_list:
        mag_name = os.path.join(data_dir, '%s_mag_grid.txt' % name)
        fit_name = os.path.join(data_dir, '%s_fit_grid.txt' % name)

        mag_data[name] = np.genfromtxt(mag_name, dtype=mag_dtype)
        fit_data[name] = np.genfromtxt(fit_name, dtype=fit_dtype)

        total_mag = mag_data[name]['ct'].sum()
        total_fit = fit_data[name]['ct'].sum()

        mag_resid_dist[name] = {}

        for cut in color_cuts:
            valid = np.where(fit_data[name]['color_resid'] <= cut)
            mag_resid_dist[name][cut] = fit_data[name][valid]


    ###plot magnitude grid
    plt.figsize = (30, 30)
    for i_fig, name in enumerate(name_list):
        plt.subplot(3, 2, i_fig+1)
        sorted_dex = np.argsort(mag_data[name]['ct'])
        plt.scatter(mag_data[name]['input'][sorted_dex],
                    mag_data[name]['mag_resid'][sorted_dex],
                    c=np.log10(mag_data[name]['ct'][sorted_dex]),
                    s=5, edgecolor='',
                    cmap=plt.cm.gist_ncar)
        plt.colorbar()
        if i_fig == 0:
            plt.xlabel('$m_{input}$')
            plt.ylabel('$m_{fit}-m_{input}$')

        xmax = mag_data[name]['input'].max()
        xmin = mag_data[name]['input'].min()
        ymax = mag_data[name]['mag_resid'].max()
        ymin = mag_data[name]['mag_resid'].min()
        plt.text(xmax-0.2*(xmax-xmin), ymax-0.1*(ymax-ymin),
                 name, fontsize=15)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'mag_validation_plot.png'))
    plt.close()


    ###plot fit grid
    plt.figsize = (30, 30)
    for i_fig, name in enumerate(name_list):
        plt.subplot(3, 2, i_fig+1)
        sorted_dex = np.argsort(fit_data[name]['ct'])
        plt.scatter(fit_data[name]['mag_resid'][sorted_dex],
                    fit_data[name]['color_resid'][sorted_dex],
                    c=np.log10(mag_data[name]['ct'][sorted_dex]),
                    s=5, edgecolor='',
                    cmap=plt.cm.gist_ncar)
        plt.colorbar()
        if i_fig == 0:
            plt.xlabel('$m_{fit} - m_{input}$')
            plt.ylabel('color residual')

        xmax = fit_data[name]['mag_resid'].max()
        xmin = fit_data[name]['mag_resid'].min()
        ymax = fit_data[name]['color_resid'].max()
        ymin = fit_data[name]['color_resid'].min()
        plt.text(xmax-0.2*(xmax-xmin), ymax-0.1*(ymax-ymin),
                 name, fontsize=15)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'color_validation_plot.png'))
    plt.close()


    ###plot distributions
    plt.figsize = (30, 30)
    text_dict = {}
    for i_fig, name in enumerate(name_list):
        plt.subplot(3, 2, i_fig+1)
        text = name

        xx, yy = get_cumulative(mag_data[name]['mag_resid'], mag_data[name]['ct'])

        five_dexes = np.where(np.abs(mag_data[name]['mag_resid'])>=5.0)
        n_five = mag_data[name]['ct'][five_dexes].sum()

        text += '\nno cut: %.1e tot, %.1e bad' % (mag_data[name]['ct'].sum(), n_five)

        xx = np.log10(xx+0.001)
        yy = np.log10(yy)

        baseline = plt.plot(xx, yy, color='k')

        xmax = xx.max()
        xmin = xx.min()
        ymax = yy.max()
        ymin = yy.min()

        cut_labels = {}

        for cut, color in zip(color_cuts, color_colors):
            xx, yy = get_cumulative(mag_resid_dist[name][cut]['mag_resid'],
                                    mag_resid_dist[name][cut]['ct'])

            five_dexes = np.where(np.abs(mag_resid_dist[name][cut]['mag_resid'])>=5.0)
            n_five = mag_resid_dist[name][cut]['ct'][five_dexes].sum()

            text += '\n$\sigma_c<=$ %.2f: %.1e tot, %.1e bad' % \
            (cut, mag_resid_dist[name][cut]['ct'].sum(), n_five)

            xx = np.log10(xx+0.001)
            yy = np.log10(yy)
            cut_labels[cut] = plt.plot(xx, yy, color=color)


        plt.text(xmax-0.95*(xmax-xmin), ymax-0.8*(ymax-ymin), text, fontsize=8)
        if i_fig == 0:
            plt.xlabel('log10($|m_{fit}-m_{input}|$)')
            plt.ylabel('log10($N_{\Delta m >=}$)')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'residual_distribution.png'))
    plt.close()
