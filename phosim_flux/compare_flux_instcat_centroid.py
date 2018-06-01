import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Sed, BandpassDict, PhotometricParameters
from lsst.sims.photUtils import getImsimFluxNorm
from lsst.sims.coordUtils import pupilCoordsFromPixelCoordsLSST
from lsst.sims.coordUtils import focalPlaneCoordsFromPupilCoordsLSST

import argparse


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

    x_mesh=np.arange(xx.min(),xx.max()+0.1,dx)
    y_mesh=np.arange(yy.min(),yy.max()+0.1,dy)
    x_mesh,y_mesh = np.meshgrid(x_mesh,y_mesh,indexing='xy')
    z_mesh = np.zeros(shape=x_mesh.shape, dtype=int)
    ct_1000b = 0

    for dex, ct in zip(dex_list, ct_list):
        ix = i_x_arr[dex]
        iy = i_y_arr[dex]
        z_mesh[iy][ix] += ct

    z_mesh = np.ma.masked_where(z_mesh==0,z_mesh)
    z_mesh = z_mesh/z_mesh.sum()
    plt.pcolormesh(x_mesh,y_mesh,z_mesh, vmin=vmin, vmax=vmax)
                   #norm=matplotlib.colors.LogNorm(vmin=1.0,
                   #                               vmax=1.2e6))
    plt.colorbar()

    ct_1000 = 0
    big_min = np.round((2.8-xx.min())/dx).astype(int)
    big_max = x_mesh.shape[0]


def process_instance_catalog(catalog_name, centroid_dir, bp_dict):

    sed_dir = getPackageDir('sims_sed_library')

    i_filter = -1
    obshistid = -1
    with open(catalog_name, 'r') as in_file:
        for line in in_file:
            params = line.strip().split()
            if params[0] == 'filter':
                i_filter = int(params[1])
            if params[0] == 'obshistid':
                obshistid = int(params[1])

            if i_filter>=0 and obshistid>=0:
                break

    filter_name = 'ugrizy'[i_filter]
    bp = bp_dict[filter_name]

    objid_arr = []
    flux_arr = []
    redshift_arr = []
    magnorm_arr = []
    phot_params = PhotometricParameters(nexp=1, exptime=30.0, gain=1.0)

    with open(catalog_name, 'r') as in_file:
        for line in in_file:
            params = line.strip().split()
            if params[0] != 'object':
                 continue

            sed_name = params[5]
            obj_id = int(params[1])
            magnorm = float(params[4])
            redshift = float(params[6])
            spec = Sed()
            spec.readSED_flambda(os.path.join(sed_dir, sed_name))
            fnorm = getImsimFluxNorm(spec, magnorm)
            spec.multiplyFluxNorm(fnorm)
            spec.redshiftSED(redshift, dimming=True)
            adu = spec.calcADU(bp, photParams=phot_params)
            objid_arr.append(obj_id)
            flux_arr.append(adu)
            redshift_arr.append(redshift)
            magnorm_arr.append(magnorm)

    centroid_files = os.listdir(centroid_dir)
    phosim_objid_arr = []
    phosim_flux_arr = []
    x_arr = []
    y_arr = []
    chip_name_arr = []
    dtype = np.dtype([('id', int), ('phot', float), ('x', float), ('y', float)])
    for file_name in centroid_files:
        if 'e_%d' % obshistid not in file_name or 'f%d' % i_filter not in file_name:
            continue
        name_params = file_name.split('_')
        raft = name_params[5]
        sensor = name_params[6]
        chip_name = '%s:%s,%s %s:%s,%s' % (raft[0],raft[1],raft[2],sensor[0],sensor[1],sensor[2])
        full_name = os.path.join(centroid_dir, file_name)
        data = np.genfromtxt(full_name, dtype=dtype, skip_header=1)
        for line in data:
            phosim_objid_arr.append(line['id'])
            phosim_flux_arr.append(line['phot'])
            x_arr.append(line['x'])
            y_arr.append(line['y'])
            chip_name_arr.append(chip_name)

    objid_arr = np.array(objid_arr)
    flux_arr = np.array(flux_arr)
    redshift_arr = np.array(redshift_arr)
    magnorm_arr = np.array(magnorm_arr)
    phosim_objid_arr = np.array(phosim_objid_arr)
    phosim_flux_arr = np.array(phosim_flux_arr)
    x_arr = np.array(x_arr)
    y_arr = np.array(y_arr)
    chip_name_arr = np.array(chip_name_arr)

    sorted_dex = np.argsort(objid_arr)
    objid_arr = objid_arr[sorted_dex]
    flux_arr = flux_arr[sorted_dex]
    redshift_arr = redshift_arr[sorted_dex]
    magnorm_arr = magnorm_arr[sorted_dex]

    sorted_dex = np.argsort(phosim_objid_arr)
    phosim_objid_arr = phosim_objid_arr[sorted_dex]
    phosim_flux_arr = phosim_flux_arr[sorted_dex]
    x_arr = x_arr[sorted_dex]
    y_arr = y_arr[sorted_dex]
    chip_name_arr = chip_name_arr[sorted_dex]

    np.testing.assert_array_equal(phosim_objid_arr, objid_arr)

    return (objid_arr, flux_arr, phosim_flux_arr,
            redshift_arr, magnorm_arr,
            x_arr, y_arr, chip_name_arr, filter_name)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--instcat_dir', type=str)
    parser.add_argument('--centroid_dir', type=str)
    parser.add_argument('--throughputs_dir', type=str)
    parser.add_argument('--fig_dir', type=str)
    args = parser.parse_args()

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles(bandpassDir=args.throughputs_dir)

    instcat_list = os.listdir(args.instcat_dir)

    phosim_flux_dict = {}
    objid_dict = {}
    catsim_flux_dict = {}
    magnorm_dict = {}
    redshift_dict = {}
    xpix_dict = {}
    ypix_dict = {}
    chip_name_dict = {}

    for instcat in instcat_list:
        instcat_file = os.path.join(args.instcat_dir, instcat)
        print('processing %s' % instcat)

        (objid, flux, phosim_flux,
         redshift, magnorm,
         xpix, ypix, chip_name, filter_name) = process_instance_catalog(instcat_file,
                                                                        args.centroid_dir,
                                                                        bp_dict)


        objid_dict[filter_name] = objid
        phosim_flux_dict[filter_name] = phosim_flux
        catsim_flux_dict[filter_name] = flux
        magnorm_dict[filter_name] = magnorm
        redshift_dict[filter_name] = redshift
        xpix_dict[filter_name] = xpix
        ypix_dict[filter_name] = ypix
        chip_name_dict[filter_name] = chip_name

        dmag = -2.5*np.log10(flux/phosim_flux)
        phosim_mag = -2.5*np.log10(phosim_flux)
        plt.figsize=(30,30)
        plt.subplot(2,2,1)
        plot_color_mesh(phosim_mag, dmag, 0.01, 0.01)
        plt.xlabel('-2.5*log10(phosim_flux)', fontsize=7)
        plt.ylabel('-2.5*log10(catsim_flux/phosim_flux)', fontsize=7)
        plt.title(filter_name, fontsize=7)

        plt.subplot(2,2,2)
        plot_color_mesh(phosim_mag, dmag, 0.01, 0.01)
        plt.xlabel('-2.5*log10(phosim_flux)', fontsize=7)
        plt.ylabel('-2.5*log10(catsim_flux/phosim_flux)', fontsize=7)
        plt.ylim(-0.5,0.5)

        plt.subplot(2,2,3)
        plot_color_mesh(redshift, dmag, 0.01, 0.01)
        plt.ylabel('-2.5*log10(catsim_flux/phosim_flux)', fontsize=7)
        plt.xlabel('redshift',fontsize=7)

        plt.subplot(2,2,4)
        plot_color_mesh(magnorm, dmag, 0.01, 0.01)
        plt.ylabel('-2.5*log10(catsim_flux/phosim_flux)', fontsize=7)
        plt.xlabel('magnorm', fontsize=7)

        plt.tight_layout()

        fig_name = os.path.join(args.fig_dir, instcat+'.png')
        plt.savefig(fig_name)
        plt.close()

    for i_f_1 in range(5):
        plt.figsize = (30,30)

        i_f_2 = i_f_1 + 1
        filter_1 = 'ugrizy'[i_f_1]
        filter_2 = 'ugrizy'[i_f_2]
        np.testing.assert_array_equal(objid_dict[filter_1], objid_dict[filter_2])
        np.testing.assert_array_equal(magnorm_dict[filter_1], magnorm_dict[filter_2])
        np.testing.assert_array_equal(redshift_dict[filter_1], redshift_dict[filter_2])

        plt.title('%s-%s' % (filter_1, filter_2), fontsize=7)
        phosim_color = 2.5*np.log10(phosim_flux_dict[filter_1]/phosim_flux_dict[filter_2])
        catsim_color = 2.5*np.log10(catsim_flux_dict[filter_1]/catsim_flux_dict[filter_2])
        catsim_color -= phosim_color
        plot_color_mesh(phosim_color, catsim_color, 0.01, 0.01)
        plt.xlabel('phosim color', fontsize=7)
        plt.ylabel('catsim_color-phosim_color', fontsize=7)

        fig_name = os.path.join(args.fig_dir,'%s_%s_color_plot.png' % (filter_1, filter_2))
        plt.savefig(fig_name)
        plt.close()

    for filter_name in 'ugrizy':
        plt.figsize = (30,30)
        plt.subplot(1,2,1)
        phosim_mag = -2.5*np.log10(phosim_flux_dict[filter_name])
        ratio = catsim_flux_dict[filter_name]/phosim_flux_dict[filter_name]
        plot_color_mesh(phosim_mag, ratio, 0.01, 0.01)
        plt.title('%s' % filter_name, fontsize=7)
        plt.xlabel('-2.5*log10(phosim_flux)')
        plt.ylabel('catsim_flux/phosim_flux')

        plt.subplot(1,2,2)
        plot_color_mesh(phosim_mag, ratio, 0.01, 0.01)
        plt.title('%s' % filter_name, fontsize=7)
        plt.xlabel('-2.5*log10(phosim_flux)')
        plt.ylabel('catsim_flux/phosim_flux')
        plt.ylim(0.8,1.2)

        plt.tight_layout()

        fig_name = os.path.join(args.fig_dir,'%s_flux_plot.png' % filter_name)
        plt.savefig(fig_name)
        plt.close()

    for i_filter in range(5):
        filter_name = 'ugrizy'[i_filter]
        filter_1 = filter_name
        filter_2 = 'ugrizy'[i_filter+1]
        xpup, ypup = pupilCoordsFromPixelCoordsLSST(xpix_dict[filter_name],
                                                    ypix_dict[filter_name],
                                                    chipName=chip_name_dict[filter_name],
                                                    band=filter_name)

        xmm, ymm = focalPlaneCoordsFromPupilCoordsLSST(xpup, ypup, band=filter_name)

        #xmm_grid = np.arange(xmm.min(),xmm.max(),0.1)
        #ymm_grid = np.arange(ymm.min(),ymm.max(),0.1)
        #mesh = np.meshgrid(xmm_grid, ymm_grid)
        #xmm_gri

        phosim_color = 2.5*np.log10(phosim_flux_dict[filter_1]/phosim_flux_dict[filter_2])
        catsim_color = 2.5*np.log10(catsim_flux_dict[filter_1]/catsim_flux_dict[filter_2])
        dcolor = catsim_color - phosim_color
        sorted_dex = np.argsort(np.abs(dcolor))
        dcolor = dcolor[sorted_dex]
        xmm = xmm[sorted_dex]
        ymm = ymm[sorted_dex]

        plt.figsize=(30,30)
        plt.scatter(xmm, ymm, c=dcolor, s=3)
        plt.title('%s-%s' % (filter_1, filter_2), fontsize=7)
        plt.colorbar(label='catsim_color-phosim_color')
        fig_name = 'focal_plane_dcolor_%s_%s.png' % (filter_1, filter_2)
        plt.savefig(os.path.join(args.fig_dir, fig_name))
        plt.close()
