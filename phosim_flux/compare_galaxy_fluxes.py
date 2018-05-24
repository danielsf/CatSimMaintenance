import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

cat_dir = 'galaxy_grid'
centroid_dir = os.path.join(cat_dir, 'output')
fig_dir = os.path.join(cat_dir, 'figs')

dt = np.dtype([('id', int), ('ra', float), ('dec', float),
               ('u', float), ('g', float), ('r', float),
               ('i', float), ('z', float), ('y', float)])

predictions = {}
for zz in np.arange(0.0, 0.8, 0.2):
    ref_name = 'galaxy_grid_ref_cat_%.1f.txt' % zz
    data = np.genfromtxt(os.path.join(cat_dir, ref_name),
                         dtype=dt)
    predictions[int(zz/0.1)] = data
    n_obj = len(data)

print(list(predictions.keys()))

phosim_data = {}
for i_filter in range(6):
    phosim_data[i_filter] = {}
    for zz in np.arange(0.0, 0.8, 0.2):
        i_zz = int(zz/0.1)
        phosim_data[i_filter][i_zz]=-1.0*np.ones(n_obj, dtype=float)

dt = np.dtype([('id', int), ('phot', float), ('x', float), ('y', float)])
centroid_files = os.listdir(centroid_dir)
for centroid_file_name in centroid_files:
    full_name = os.path.join(centroid_dir, centroid_file_name)
    data = np.genfromtxt(full_name, dtype=dt, skip_header=1)
    params = centroid_file_name.split('_')
    obshistid = int(params[3])
    i_filter = obshistid//100
    i_zz = obshistid-100*i_filter-1

    instcat_name = 'instcat_galaxy_grid_%.1f_%d.txt' % (i_zz*0.1, i_filter)
    valid_icat = False
    with open(os.path.join(cat_dir,instcat_name), 'r') as in_file:
        for line in in_file:
            i_params = line.strip().split()
            if i_params[0] == 'obshistid':
                assert int(i_params[1]) == obshistid
                valid_icat = True
                break
    assert valid_icat

    phosim_data[i_filter][i_zz][data['id']-1] = data['phot']


for i_filter_1 in range(5):
    i_filter_2 = i_filter_1+1
    filter_1 = 'ugrizy'[i_filter_1]
    filter_2 = 'ugrizy'[i_filter_2]
    true_colors = []
    delta_colors = []
    zz_arr = []

    for i_zz in range(0,8,2):
        valid = np.where(np.logical_and(phosim_data[i_filter_1][i_zz]>0.0,
                                        phosim_data[i_filter_2][i_zz]>0.0))
        for dex in valid[0]:
            c_color = predictions[i_zz][filter_1][dex]/predictions[i_zz][filter_2][dex]
            c_color = -2.5*np.log10(c_color)
            true_colors.append(c_color)
            p_color = phosim_data[i_filter_1][i_zz][dex]/phosim_data[i_filter_2][i_zz][dex]
            p_color = -2.5*np.log10(p_color)
            delta_colors.append(p_color-c_color)
            zz_arr.append(i_zz)

    true_colors = np.array(true_colors)
    delta_colors = np.array(delta_colors)
    zz_arr = np.array(zz_arr)

    plt.figsize=(30,30)
    plt.subplot(2,2,1)
    plt.title('all z', fontsize=7)
    plt.scatter(true_colors, delta_colors, s=5)
    plt.xlabel('%s-%s (CatSim)' % (filter_1, filter_2))
    plt.ylabel('PhoSim-CatSim')
    plt.axhline(0.0, linestyle='--', color='r')
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)

    i_fig_list = [2, 4, 5, 6]

    for i_fig, i_zz in zip(i_fig_list, range(0,8,2)):
        plt.subplot(3,2,i_fig)
        valid = np.where(zz_arr==i_zz)
        if len(valid[0])==0:
            continue

        plt.scatter(true_colors[valid], delta_colors[valid], s=5)

        x_min = true_colors[valid].min()
        x_max = true_colors[valid].max()
        y_min = delta_colors[valid].min()
        y_max = delta_colors[valid].max()
        if 0.0<y_min:
            y_min=0.0
        if 0.0>y_max:
            y_max=0.0

        dx=x_max-x_min
        dy=y_max-y_min
        corner_x = np.array([x_min+0.25*dx,x_min+0.25*dx, x_max-0.25*dx,x_max-0.25*dx])
        corner_y = np.array([y_min+0.25*dy,y_max-0.25*dy,y_min+0.25*dy,y_max-0.25*dy])
        i_best = -1
        d_best = 0.0
        for i_c in range(4):
            dd = np.sqrt((true_colors[valid]-corner_x[i_c])**2+
                         (delta_colors[valid]-corner_y[i_c])**2)

            dd_min = np.min(dd)
            if dd_min>d_best:
                d_best=dd_min
                i_best=i_c

        plt.text(corner_x[i_best],corner_y[i_best],
                 'z=%.1f' % (0.1*i_zz),fontsize=7)

        #plt.xlabel('%s-%s (CatSim)' % (filter_1, filter_2))
        #plt.ylabel('PhoSim-CatSim')
        plt.axhline(0.0, linestyle='--', color='r')
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)

    fig_name = os.path.join(fig_dir,'colors_%s_%s.eps' % (filter_1, filter_2))
    #plt.tight_layout()
    plt.savefig(fig_name)
    plt.close()

for i_filter in range(6):
    filter_name = 'ugrizy'[i_filter]
    true_fluxes = []
    flux_ratios = []
    for i_zz in range(0,8,2):
        valid = np.where(phosim_data[i_filter][i_zz]>0.0)
        for dex in valid[0]:
            c_flux = predictions[i_zz][filter_name][dex]
            p_flux = phosim_data[i_filter][i_zz][dex]
            true_fluxes.append(c_flux)
            flux_ratios.append(p_flux/c_flux)

    true_fluxes = np.array(true_fluxes)
    flux_ratios = np.array(flux_ratios)

    plt.figsize = (30,30)
    plt.scatter(true_fluxes, flux_ratios, s=5)
    plt.xlabel('%s ADU (CatSim)' % filter_name)
    plt.ylabel('PhoSim/CatSim (should be == gain)')
    plt.axhline(2.3,linestyle='--', color='g')
    plt.axhline(1.0, linestyle='--', color='r')
    plt.yscale('log')
    plt.xscale('log')
    fig_name = os.path.join(fig_dir, 'fluxes_%s.eps' % filter_name)
    plt.savefig(fig_name)
    plt.close()
