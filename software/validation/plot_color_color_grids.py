from __future__ import with_statement
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

if __name__ == "__main__":

    fig_name_list = [os.path.join(args.output_dir, '%s_fit_grid.png' % args.prefix),
                     os.path.join(args.output_dir, '%s_input_grid.png' % args.prefix)]
    grid_list = [fit_grids, input_grids]

    for fig_name, grid_set in zip(fig_name_list, grid_list):
        plt.figsize = (30, 30)
        for i_fig in range(len(mags)-2):
            plt.subplot(2, 2, i_fig+1)
            mag1 = mags[i_fig]
            mag2 = mags[i_fig+1]
            mag3 = mags[i_fig+2]

            grid =grid_set[mag1]

            x_arr = []
            y_arr = []
            ct_arr = []
            for xx in grid:
                for yy in grid[xx]:
                    x_arr.append(xx*dmag)
                    y_arr.append(yy*dmag)
                    ct_arr.append(grid[xx][yy])
            x_arr = np.array(x_arr)
            y_arr = np.array(y_arr)
            ct_arr = np.array(ct_arr)

            if len(ct_arr)==0:
                continue

            sorted_dexes = np.argsort(ct_arr)
            plt.scatter(x_arr[sorted_dexes], y_arr[sorted_dexes],
                        c=ct_arr[sorted_dexes],
                        edgecolor='',s=5,
                        cmap=plt.cm.gist_ncar,
                        norm=Normalize(vmin=ct_arr.min(), vmax=ct_arr.max()))

            cb = plt.colorbar()

            plt.xlabel('%s-%s' % (mag1, mag2))
            plt.ylabel('%s-%s' % (mag2, mag3))

        plt.tight_layout()
        plt.savefig(fig_name)
        plt.close()
