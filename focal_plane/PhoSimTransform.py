import numpy as np
import os
import palpy

from lsst.sims.utils.CodeUtilities import _validate_inputs

__all__ = ["PhoSimPixelTransformer",
           "_rawPupilCoordsFromObserved",
           "_rawObservedFromPupilCoords"]


def load_focal_plane(perturbed=True):
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
            if perturbed:
                local_dict['rot'] = float(params[10])
            else:
                local_dict['rot'] = 0.0
            assert np.abs(float(params[11]))<1.0e-10
            assert np.abs(float(params[12]))<1.0e-10
            if perturbed:
                local_dict['dx'] = float(params[13])
                local_dict['dy'] = float(params[14])
            else:
                local_dict['dx'] = 0.0
                local_dict['dy'] = 0.0
            assert np.abs(float(params[15]))<1.0e-10
            chip_data[name] = local_dict

    return chip_data


class PhoSimPixelTransformer(object):

    def __init__(self, perturbed=True):
        self._chip_data = load_focal_plane(perturbed=perturbed)

    def _chip_center_mm(self, chipName):
        chip= self._chip_data[chipName]
        x0 = 0.001*chip['x']+chip['dx']
        y0 = 0.001*chip['y']+chip['dy']
        return x0, y0

    def _chip_rot_matrix(self, chipName):
        chip = self._chip_data[chipName]

        theta = np.radians(chip['rot'])
        cc = np.cos(theta)
        ss = np.sin(theta)
        rotMat = np.array([[cc, -ss], [ss, cc]])
        return rotMat


    def mmFromPix(self, xpix, ypix, chipName):
        chip = self._chip_data[chipName]
        x0, y0 = self._chip_center_mm(chipName)
        dp = 0.001*chip['p_size']

        # rotate (1,0) and (0,1) so that x_vec and y_vec
        # reflect actual orientation of the pixel grid
        rotMat = self._chip_rot_matrix(chipName)

        x_pix_vec = np.array([1.0, 0.0])
        y_pix_vec = np.array([0.0, 1.0])
        x_pix_vec = np.dot(rotMat, x_pix_vec)
        y_pix_vec = np.dot(rotMat, y_pix_vec)

        mm_vec = np.array([x0, y0])
        mm_vec += (xpix+0.5-chip['n_x']/2)*dp*x_pix_vec
        mm_vec += (ypix+0.5-chip['n_y']/2)*dp*y_pix_vec

        return mm_vec[0], mm_vec[1]

    def pixFromMM(self, xmm, ymm, chipName):
        """
        xmm -- x coord in millimeters
        ymm -- y coord in millimeters
        chipName = like R:2,2 S:1,1
        """
        x0, y0 = self._chip_center_mm(chipName)

        chip = self._chip_data[chipName]
        dp = 0.001*chip['p_size']
        xpix0 = -0.5+chip['n_x']/2
        ypix0 = -0.5+chip['n_y']/2

        rotMat = self._chip_rot_matrix(chipName)

        x_pix_vec = np.array([1.0, 0.0])
        y_pix_vec = np.array([0.0, 1.0])
        x_pix_vec = np.dot(rotMat, x_pix_vec)
        y_pix_vec = np.dot(rotMat, y_pix_vec)

        d_mm_vec = np.array([xmm-x0, ymm-y0])

        dx = np.dot(x_pix_vec, d_mm_vec)
        dy = np.dot(y_pix_vec, d_mm_vec)

        xpix = xpix0 + dx/dp
        ypix = ypix0 + dy/dp

        return xpix, ypix


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
