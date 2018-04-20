import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

def make_2d_histogram(xx, yy, dx, dy):
    """
    returns indices and counts of unique points on the map
    """
    i_color1 = np.round(xx/dx).astype(int)
    i_color2 = np.round(yy/dy).astype(int)
    dex_reverse = np.array([i_color1, i_color2])
    dex_arr = dex_reverse.transpose()
    # see http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
    dex_raw = np.ascontiguousarray(dex_arr).view(np.dtype((np.void, dex_arr.dtype.itemsize*dex_arr.shape[1])))
    _, unique_rows, unique_counts = np.unique(dex_raw, return_index=True, return_counts=True)

    return unique_rows, unique_counts


def plot_color(xx, yy, dx, dy):
    dexes, cts = make_2d_histogram(xx, yy, dx, dy)
    sorted_dex = np.argsort(cts)
    dexes = dexes[sorted_dex]
    cts = cts[sorted_dex]
    plt.scatter(xx[dexes], yy[dexes], c=cts, s=5,
                cmap=plt.cm.gist_ncar, edgecolor='')

    plt.colorbar()


def plot_color_mesh(xx, yy, dx, dy, vmin=None, vmax=None):
    i_x_arr = np.round((xx-xx.min())/dx).astype(int)
    i_y_arr = np.round((yy-yy.min())/dy).astype(int)
    new_x = i_x_arr*dx
    new_y = i_y_arr*dy
    dex_list, ct_list = make_2d_histogram(new_x, new_y, dx, dy)

    if i_x_arr.min()<0 or i_y_arr.min()<0:
        raise RuntimeError('negative dex')

    x_mesh=np.arange(xx.min(),xx.max()+dx,dx)
    y_mesh=np.arange(yy.min(),yy.max()+dy,dy)
    x_mesh,y_mesh = np.meshgrid(x_mesh,y_mesh,indexing='xy')
    z_mesh = np.zeros(shape=x_mesh.shape, dtype=int)
    ct_1000b = 0

    for dex, ct in zip(dex_list, ct_list):
        ix = i_x_arr[dex]
        iy = i_y_arr[dex]
        z_mesh[iy][ix] += ct

    z_mesh = np.ma.masked_where(z_mesh==0,z_mesh)
    plt.pcolormesh(x_mesh,y_mesh,z_mesh, vmin=vmin, vmax=vmax)
                   #norm=matplotlib.colors.LogNorm(vmin=1.0,
                   #                               vmax=1.2e6))

    ct_1000 = 0
    big_min = np.round((2.8-xx.min())/dx).astype(int)
    big_max = x_mesh.shape[0]


if __name__ == "__main__":

    dtype = np.dtype([('redshift', float), ('u', float),
                      ('g', float), ('r', float),
                      ('i', float), ('z', float),
                      ('y', float),
                      ('tau', float), ('sfu', float),
                      ('sfg', float), ('sfr', float),
                      ('sfi', float), ('sfz', float),
                      ('sfy', float)])

    data = np.genfromtxt('agn_variability_distribution.txt',
                         dtype=dtype)

    plt.figsize = (30,30)
    for i_band, band in enumerate('ugrizy'):
        plt.subplot(3,2,i_band+1)

        valid = np.where(data['redshift']<=4.0)

        plot_color_mesh(data['redshift'][valid], data[band][valid], 0.05, 0.1)
        if i_band==0:
            plt.xlabel('redshift', fontsize=7)

        plt.ylabel('%s (obs)' % band, fontsize=7)
        if i_band==0:
            plt.colorbar(label='sources per pixel')

        yticks = np.arange(np.floor(data[band].min()), np.ceil(data[band].max()),
                           2)

        ylabels = ['' if ii%2==0 else '%d' % yticks[ii]
                   for ii in range(len(yticks))]

        plt.yticks(yticks, ylabels)
    plt.tight_layout()
    plt.savefig('quiescent_mag_vs_redshift.png')
    plt.close()

    plt.figsize = (30,30)
    for i_band, band in enumerate('ugrizy'):
        sf_name = 'sf%s' % band
        plt.subplot(3,2,i_band+1)

        valid = np.where(data['redshift']<=4.0)

        plot_color_mesh(data['redshift'][valid], data[sf_name][valid], 0.05, 0.05)
        if i_band==0:
            plt.xlabel('redshift', fontsize=7)

        plt.ylabel('SF(%s)' % band, fontsize=7)
        if i_band==0:
            plt.colorbar(label='sources per pixel')

        yticks = np.arange(np.floor(data[sf_name][valid].min()),
                           np.ceil(data[sf_name][valid].max()),
                           0.25)

        ylabels = ['' if ii%3!=0 else '%.2f' % yticks[ii]
                   for ii in range(len(yticks))]

        plt.yticks(yticks, ylabels)
    plt.tight_layout()
    plt.savefig('structure_function_vs_redshift.png')
    plt.close()

    plt.figsize = (30,30)
    plt.suptitle('redshift <= 4.0')
    for i_band, band in enumerate('ugrizy'):
        sf_name = 'sf%s' % band
        plt.subplot(3,2,i_band+1)

        valid = np.where(data['redshift']<=4.0)

        plot_color_mesh(data[band][valid], data[sf_name][valid], 0.05, 0.05)
        plt.xlabel('%s (obs)' % band, fontsize=7)

        plt.ylabel('SF(%s)' % band, fontsize=7)
        if i_band==0:
            plt.colorbar(label='sources per pixel')

        yticks = np.arange(np.floor(data[sf_name][valid].min()),
                           np.ceil(data[sf_name][valid].max()),
                           0.25)

        ylabels = ['' if ii%3!=0 else '%.2f' % yticks[ii]
                   for ii in range(len(yticks))]

        plt.yticks(yticks, ylabels)
    plt.tight_layout(rect=[0,0,1,0.95])
    plt.savefig('structure_function_vs_quiescent_magnitude.png')
    plt.close()

    plt.figsize = (30,30)
    plt.suptitle('redshift <= 4.0')
    for i_band, band in enumerate('ugrizy'):
        sf_name = 'sf%s' % band
        plt.subplot(3,2,i_band+1)

        valid = np.where(data['redshift']<=4.0)

        tau_obs = np.log10(data['tau'][valid]/(1.0+data['redshift'][valid]))

        plot_color_mesh(data[sf_name][valid],
                        tau_obs,
                        0.05, 0.1)

        plt.xlabel('SF(%s)' % band, fontsize=7)
        if i_band==0:
            plt.colorbar(label='sources per pixel')
            plt.ylabel('log10[tau (days)]', fontsize=7)

        xticks = np.arange(np.floor(data[sf_name][valid].min()),
                           np.ceil(data[sf_name][valid].max()),
                           0.25)

        xlabels = ['' if ii%3!=0 else '%.2f' % xticks[ii]
                   for ii in range(len(xticks))]

        plt.xticks(xticks, xlabels)

        yticks = np.arange(0.0, 4.0, 0.25)
        ylabels = ['' if ii%4!=0 else '%d' % yticks[ii]
                   for ii in range(len(yticks))]
        plt.yticks(yticks, ylabels)

    plt.tight_layout(rect=[0,0,1,0.95])
    plt.savefig('structure_function_vs_tau.png')
    plt.close()

    sf_name = None

    plt.figsize = (30,30)
    plt.suptitle('redshift <= 4.0')
    for i_band, band in enumerate('ugrizy'):
        plt.subplot(3,2,i_band+1)

        valid = np.where(data['redshift']<=4.0)

        tau_obs = np.log10(data['tau'][valid]/(1.0+data['redshift'][valid]))

        plot_color_mesh(data[band][valid],
                        tau_obs,
                        0.05, 0.1)

        plt.xlabel('%s (obs)' % band, fontsize=7)
        if i_band==0:
            plt.colorbar(label='sources per pixel')
            plt.ylabel('log10[tau (days)]', fontsize=7)

        xticks = np.arange(np.floor(data[band][valid].min()),
                           np.ceil(data[band][valid].max()),
                           1.0)

        xlabels = ['' if ii%5!=0 else '%d' % xticks[ii]
                   for ii in range(len(xticks))]

        plt.xticks(xticks, xlabels)


        yticks = np.arange(0.0, 4.0, 0.25)
        ylabels = ['' if ii%4!=0 else '%d' % yticks[ii]
                   for ii in range(len(yticks))]
        plt.yticks(yticks, ylabels)

    plt.tight_layout(rect=[0,0,1,0.95])
    plt.savefig('tau_vs_quiescent_magnitude.png')
    plt.close()
