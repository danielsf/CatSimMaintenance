import os
import copy
import numpy as np
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap
from lsst.sims.catUtils.exampleCatalogDefinitions import write_phoSim_header
from lsst.sims.catUtils.mixins import PhoSimAstrometryStars
from lsst.sims.utils import _observedFromICRS
from lsst.sims.coordUtils import lsst_camera
from lsst.sims.coordUtils import raDecFromPixelCoords
from lsst.sims.coordUtils import pixelCoordsFromRaDecLSST
from lsst.sims.coordUtils import focalPlaneCoordsFromPupilCoords
from lsst.sims.coordUtils import pupilCoordsFromFocalPlaneCoords
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST
from lsst.sims.utils import Site
from lsst.afw.cameraGeom import SCIENCE
from lsst.sims.GalSimInterface import LSSTCameraWrapper
from lsst.sims.utils.CodeUtilities import _validate_inputs
from lsst.sims.utils import arcsecFromRadians
import palpy

import argparse

def _rawPupilCoordsFromObserved(ra_obs, dec_obs, ra0, dec0, rotSkyPos):
    """
    Convert Observed RA, Dec into pupil coordinates

    Parameters
    ----------
    ra_obs is the observed RA in radians

    dec_obs is the observed Dec in radians

    obs_metadata is an ObservationMetaData characterizing the telescope location and pointing

    epoch is the epoch of the mean RA and Dec in julian years (default=2000.0)

    includeRefraction is a boolean controlling the application of refraction.

    Returns
    --------
    A numpy array whose first row is the x coordinate on the pupil in
    radians and whose second row is the y coordinate in radians
    """

    are_arrays = _validate_inputs([ra_obs, dec_obs], ['ra_obs', 'dec_obs'],
                                  "pupilCoordsFromObserved")

    theta = -1.0*rotSkyPos

    ra_pointing = ra0
    dec_pointing = dec0

    # palpy.ds2tp performs the gnomonic projection on ra_in and dec_in
    # with a tangent point at (pointingRA, pointingDec)
    #
    if not are_arrays:
        try:
            x, y = palpy.ds2tp(ra_obs, dec_obs, ra_pointing, dec_pointing)
        except:
            x = np.NaN
            y = np.NaN
    else:
        try:
            x, y = palpy.ds2tpVector(ra_obs, dec_obs, ra_pointing, dec_pointing)
        except:
            # apparently, one of your ra/dec values was improper; we will have to do this
            # element-wise, putting NaN in the place of the bad values
            x = []
            y = []
            for rr, dd in zip(ra_obs, dec_obs):
                try:
                    xx, yy = palpy.ds2tp(rr, dd, ra_pointing, dec_pointing)
                except:
                    xx = np.NaN
                    yy = np.NaN
                x.append(xx)
                y.append(yy)
            x = np.array(x)
            y = np.array(y)

    # rotate the result by rotskypos (rotskypos being "the angle of the sky relative to
    # camera coordinates" according to phoSim documentation) to account for
    # the rotation of the focal plane about the telescope pointing

    x_out = x*np.cos(theta) - y*np.sin(theta)
    y_out = x*np.sin(theta) + y*np.cos(theta)

    return np.array([x_out, y_out])

def _rawObservedFromPupilCoords(xPupil, yPupil, ra0, dec0, rotSkyPos):
    """
    Convert pupil coordinates into observed (RA, Dec)

    @param [in] xPupil -- pupil coordinates in radians.
    Can be a numpy array or a number.

    @param [in] yPupil -- pupil coordinates in radians.
    Can be a numpy array or a number.

    @param [in] obs_metadata -- an instantiation of ObservationMetaData characterizing
    the state of the telescope

    @param [in] epoch -- julian epoch of the mean equinox used for the coordinate
    transformations (in years; defaults to 2000)

    @param[in] includeRefraction -- a boolean which controls the effects of refraction
    (refraction is used when finding the observed coordinates of the boresite specified
    by obs_metadata)

    @param [out] a 2-D numpy array in which the first row is observed RA and the second
    row is observed Dec (both in radians).  Note: these are not ICRS coordinates.
    These are RA and Dec-like coordinates resulting from applying precession, nutation,
    diurnal aberration and annual aberration on top of ICRS coordinates.

    WARNING: This method does not account for apparent motion due to parallax.
    This method is only useful for mapping positions on a theoretical focal plane
    to positions on the celestial sphere.
    """

    are_arrays = _validate_inputs([xPupil, yPupil], ['xPupil', 'yPupil'],
                                  "observedFromPupilCoords")

    ra_pointing = ra0
    dec_pointing = dec0

    # This is the same as theta in pupilCoordsFromRaDec, except without the minus sign.
    # This is because we will be reversing the rotation performed in that other method.
    theta = rotSkyPos

    x_g = xPupil*np.cos(theta) - yPupil*np.sin(theta)
    y_g = xPupil*np.sin(theta) + yPupil*np.cos(theta)

    # x_g and y_g are now the x and y coordinates
    # can now use the PALPY method palDtp2s to convert to RA, Dec.

    if are_arrays:
        raObs, decObs = palpy.dtp2sVector(x_g, y_g, ra_pointing, dec_pointing)
    else:
        raObs, decObs = palpy.dtp2s(x_g, y_g, ra_pointing, dec_pointing)

    return raObs, decObs


def load_focal_plane():
    phosim_dir = os.path.join("/Users/danielsf/physics/phosim_release")
    phosim_file = os.path.join(phosim_dir, "data", "lsst", "focalplanelayout.txt")

    chip_data = {}
    with open(phosim_file, 'r') as input_file:
        for line in input_file:
            if line[0] == '#':
                continue
            params = line.strip().split()
            mangled_name = params[0]
            name = mangled_name[0]+':'+mangled_name[1]+','+mangled_name[2]+' '
            name += mangled_name[4]+':'+mangled_name[5]+','+mangled_name[6]
            local_dict = {}
            local_dict['x'] = float(params[1])
            local_dict['y'] = float(params[2])
            local_dict['p_size'] = float(params[3])
            local_dict['n_x'] = int(params[4])
            local_dict['n_y'] = int(params[5])
            local_dict['rot'] = float(params[10])
            assert np.abs(float(params[11]))<1.0e-10
            assert np.abs(float(params[12]))<1.0e-10
            local_dict['dx'] = float(params[13])
            local_dict['dy'] = float(params[14])
            assert np.abs(float(params[15]))<1.0e-10
            chip_data[name] = local_dict

    return chip_data


class PhoSimPixelTransformer(object):

    def __init__(self):
        self._chip_data = load_focal_plane()

    def mmFromPix(self, xpix, ypix, chipName):
        chip = self._chip_data[chipName]
        x0 = 0.001*chip['x']+chip['dx']
        y0 = 0.001*chip['y']+chip['dy']
        dp = 0.001*chip['p_size']

        # rotate (1,0) and (0,1) so that x_vec and y_vec
        # reflect actual orientation of the pixel grid
        theta = np.radians(chip['rot'])
        cc = np.cos(theta)
        ss = np.sin(theta)
        rotMat = np.array([[cc, -ss], [ss, cc]])

        x_vec = np.array([1.0, 0.0])
        y_vec = np.array([0.0, 1.0])
        x_vec = np.dot(rotMat, x_vec)
        y_vec = np.dot(rotMat, y_vec)

        mm_vec = np.array([x0, y0])
        mm_vec += (xpix-chip['n_x']/2)*dp*x_vec
        mm_vec += (ypix-chip['n_y']/2)*dp*y_vec

        return mm_vec[0], mm_vec[1]


if __name__ == "__main__":

    coord_converter = PhoSimPixelTransformer()

    parser = argparse.ArgumentParser()
    parser.add_argument('--obs', type=int, default=230)
    args = parser.parse_args()

    chip_name = 'R:2,2 S:1,1'
    mangled_name = chip_name.replace(':','').replace(',','').replace(' ','_')

    opsimdb = os.path.join('/Users', 'danielsf', 'physics', 'lsst_150412',
                           'Development', 'garage', 'OpSimData',
                           'minion_1016_sqlite.db')

    assert os.path.exists(opsimdb)
    obs_gen = ObservationMetaDataGenerator(database=opsimdb)
    obs_list = obs_gen.getObservationMetaData(obsHistID=args.obs)
    obs = obs_list[0]
    site_no_atm = Site(name="LSST",
                       pressure=0.0,
                       humidity=0.0)
    obs.site=site_no_atm
    assert np.abs(obs.site.pressure)<1.0e-6
    assert np.abs(obs.site.humidity)<1.0e-6

    camera = lsst_camera()
    x_pix_arr = np.arange(100.0, 3900.0, 200.0)
    y_pix_arr = np.arange(100.0, 3900.0, 200.0)
    pix_grid = np.meshgrid(x_pix_arr, y_pix_arr)
    x_pix_arr = pix_grid[0].flatten()
    y_pix_arr = pix_grid[1].flatten()

    id_grid = np.arange(len(x_pix_arr)+1).astype(int)
    x_mm = []
    y_mm = []
    for xx, yy in zip(x_pix_arr, y_pix_arr):
        xm, ym = coord_converter.mmFromPix(xx, yy, chip_name)
        x_mm.append(xm)
        y_mm.append(ym)

    chip_data = coord_converter._chip_data[chip_name]
    x_mm.append(chip_data['x']*0.001+chip_data['dx'])
    y_mm.append(chip_data['y']*0.001+chip_data['dy'])
    x_mm = np.array(x_mm)
    y_mm = np.array(y_mm)

    x_pup_arr, y_pup_arr = pupilCoordsFromFocalPlaneCoords(x_mm, y_mm, camera=lsst_camera())
    
    camera_wrapper = LSSTCameraWrapper()

    ra_rad, dec_rad = _rawObservedFromPupilCoords(x_pup_arr, y_pup_arr,
                                                  obs._pointingRA,
                                                  obs._pointingDec,
                                                  obs._rotSkyPos)

    id_grid = np.arange(len(ra_rad))

    ra_grid = np.degrees(ra_rad)
    dec_grid = np.degrees(dec_rad)

    phosim_header_map = copy.deepcopy(DefaultPhoSimHeaderMap)
    phosim_header_map['rawSeeing'] = ('rawSeeing', None)
    phosim_header_map['FWHMeff'] = ('FWHMeff', None)
    phosim_header_map['FWHMgeom'] = ('FWHMgeom',None)
    phosim_header_map['nsnap'] = 1
    phosim_header_map['vistime'] =30.0

    with open('raw_catalogs/star_cross_%d_%s.txt' % (args.obs, mangled_name), 'w') as out_file:
        write_phoSim_header(obs, out_file, phosim_header_map)
        for (i_obj, ra, dec) in zip(id_grid, ra_grid, dec_grid):
            out_file.write('object %d ' % (i_obj))
            out_file.write('%.17f %.17f ' % (ra, dec))
            out_file.write('21.0 starSED/kurucz/km10_5750.fits_g10_5750.gz 0 0 0 0 0 0 point none CCM 0.03380581 3.1\n')

        #out_file.write('object 25 ')
        #out_file.write('%.17f %.17f ' % (obs.pointingRA, obs.pointingDec))
        #out_file.write('19.0 starSED/kurucz/km10_5750.fits_g10_5750.gz 0 0 0 0 0 0 point none CCM 0.03380581 3.1\n')

    with open('raw_catalogs/star_cross_predicted_%d_%s.txt' % (args.obs, mangled_name), 'w') as out_file:
        out_file.write('# id xpup ypup\n')
        for ii in range(len(x_pix_arr)):
            out_file.write('%d %.17f %.17f %e %e\n' %
                           (id_grid[ii],
                            x_pix_arr[ii],
                            y_pix_arr[ii],
                            x_mm[ii], y_mm[ii]))
        #out_file.write('25 R22S11 %.17f %.17f\n' % (x_pix_c, y_pix_c))
