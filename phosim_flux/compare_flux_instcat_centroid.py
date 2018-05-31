import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Sed, BandpassDict, PhotometricParameters
from lsst.sims.photUtils import getImsimFluxNorm

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
    dtype = np.dtype([('id', int), ('phot', float), ('x', float), ('y', float)])
    for file_name in centroid_files:
        if 'e_%d' % obshistid not in file_name or 'f%d' % i_filter not in file_name:
            continue
        full_name = os.path.join(centroid_dir, file_name)
        data = np.genfromtxt(full_name, dtype=dtype, skip_header=1)
        for line in data:
            phosim_objid_arr.append(line['id'])
            phosim_flux_arr.append(line['phot'])

    objid_arr = np.array(objid_arr)
    flux_arr = np.array(flux_arr)
    redshift_arr = np.array(redshift_arr)
    magnorm_arr = np.array(magnorm_arr)
    phosim_objid_arr = np.array(phosim_objid_arr)
    phosim_flux_arr = np.array(phosim_flux_arr)

    sorted_dex = np.argsort(objid_arr)
    objid_arr = objid_arr[sorted_dex]
    flux_arr = flux_arr[sorted_dex]
    redshift_arr = redshift_arr[sorted_dex]
    magnorm_arr = magnorm_arr[sorted_dex]

    sorted_dex = np.argsort(phosim_objid_arr)
    phosim_objid_arr = phosim_objid_arr[sorted_dex]
    phosim_flux_arr = phosim_flux_arr[sorted_dex]

    np.testing.assert_array_equal(phosim_objid_arr, objid_arr)

    return (objid_arr, flux_arr, phosim_flux_arr, redshift_arr, magnorm_arr)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--instcat_dir', type=str)
    parser.add_argument('--centroid_dir', type=str)
    parser.add_argument('--throughputs_dir', type=str)
    parser.add_argument('--fig_dir', type=str)
    args = parser.parse_args()

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles(bandpassDir=args.throughputs_dir)

    instcat_list = os.listdir(args.instcat_dir)

    for instcat in instcat_list:
        instcat_file = os.path.join(args.instcat_dir, instcat)
        print('processing %s' % instcat)

        (objid, flux, phosim_flux,
         redshift, magnorm) = process_instance_catalog(instcat_file, args.centroid_dir,
                                                       bp_dict)

        dmag = -2.5*np.log10(flux/phosim_flux)
        phosim_mag = -2.5*np.log10(phosim_flux)
        plt.figsize=(30,30)
        plt.subplot(2,2,1)
        plot_color_mesh(phosim_mag, dmag, 0.01, 0.01)
        plt.xlabel('-2.5*log10(phosim_flux)', fontsize=7)
        plt.ylabel('-2.5*log10(catsim_flux/phosim_flux)', fontsize=7)

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
