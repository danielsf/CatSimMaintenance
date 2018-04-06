"""
Generate an InstanceCatalog using all of the CatSim-PhoSim infrastructure.
Predict pixel positions using that same infrastructure.
"""

import os
import numpy as np
import copy
import argparse

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogPoint
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap

from lsst.sims.coordUtils import raDecFromPixelCoords
from lsst.sims.catUtils.mixins import PhoSimAstrometryBase
from lsst.sims.coordUtils import lsst_camera
from lsst.sims.catUtils.exampleCatalogDefinitions import write_phoSim_header
from lsst.sims.utils import observedFromICRS
from lsst.sims.utils import Site
from lsst.sims.GalSimInterface import LSSTCameraWrapper
from lsst.sims.coordUtils import focalPlaneCoordsFromRaDec

from lsst.afw.cameraGeom import SCIENCE

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--obs', type=int, default=230)
    args = parser.parse_args()

    camera = lsst_camera()
    det_name_list = []
    for det in camera:
        if det.getType() != SCIENCE:
            continue
        det_name_list.append(det.getName())
    det_name_list.sort()

    opsimdb = os.path.join('/Users', 'danielsf', 'physics', 'lsst_150412',
                           'Development', 'garage', 'OpSimData',
                           'minion_1016_sqlite.db')

    assert os.path.exists(opsimdb)
    obs_gen = ObservationMetaDataGenerator(database=opsimdb)
    obs_list = obs_gen.getObservationMetaData(obsHistID=args.obs)
    obs = obs_list[0]
    filter_name = obs.bandpass
    site_no_atm = Site(name="LSST",
                       pressure=0.0,
                       humidity=0.0)
    obs.site=site_no_atm
    assert np.abs(obs.site.pressure)<1.0e-6
    assert np.abs(obs.site.humidity)<1.0e-6

    xpix_0 = np.arange(200.0, 3800.0, 1000.0)
    ypix_0 = np.arange(200.0, 3800.0, 1000.0)
    pix_grid = np.meshgrid(xpix_0, ypix_0)
    cam_xpix_in = pix_grid[0].flatten()
    cam_ypix_in = pix_grid[1].flatten()

    camera_wrapper = LSSTCameraWrapper()

    de_precessor = PhoSimAstrometryBase()
    phosim_header_map = copy.deepcopy(DefaultPhoSimHeaderMap)
    phosim_header_map['nsnap'] = 1
    phosim_header_map['exptime'] = 30.0

    i_obj = 0

    with open('full_catalogs/star_grid_%d_%s.txt' % (args.obs, filter_name), 'w') as cat_file:
        write_phoSim_header(obs, cat_file, phosim_header_map)
        with open('full_catalogs/star_predicted_%d_%s.txt' % (args.obs, filter_name), 'w') as truth_file:
            truth_file.write('# i_obj ra_icrs dec_icrs ra_deprecessed dec_deprecessed x_dm y_dm x_f y_f x_cam y_cam\n')

            for det_name in det_name_list:
                det_name_m = det_name.replace(':','').replace(',','').replace(' ','_')


                dm_xpix_in, dm_ypix_in = camera_wrapper.dmPixFromCameraPix(cam_xpix_in,
                                                                           cam_ypix_in,
                                                                           det_name)

                ra_icrs, dec_icrs = raDecFromPixelCoords(dm_xpix_in, dm_ypix_in, det_name,
                                                         camera=lsst_camera(),
                                                         obs_metadata=obs)

                ra_obs, dec_obs = observedFromICRS(ra_icrs, dec_icrs, obs_metadata=obs,
                                                   epoch=2000.0, includeRefraction=False)

                (ra_deprecessed,
                 dec_deprecessed) = de_precessor._dePrecess(np.radians(ra_obs),
                                                            np.radians(dec_obs),
                                                            obs)

                x_f, y_f = focalPlaneCoordsFromRaDec(ra_icrs, dec_icrs, obs_metadata=obs,
                                                     camera=lsst_camera())

                for ra, dec, camxx, camyy, r_icrs, d_icrs, dmxx, dmyy, xxf, yyf in \
                zip(ra_deprecessed, dec_deprecessed, cam_xpix_in, cam_ypix_in, ra_icrs, dec_icrs,
                    dm_xpix_in, dm_ypix_in, x_f, y_f):

                    i_obj += 1
                    cat_file.write('object %d ' % (i_obj))
                    cat_file.write('%.17f %.17f ' % (np.degrees(ra), np.degrees(dec)))
                    cat_file.write('21.0 flatSED/sed_flat_short.txt.gz 0 0 0 0 0 0 point none CCM 0.03380581 3.1\n')

                    truth_file.write('%d %.10f %.10f %.10f %.10f %.2f %.2f %.6f %.6f %.5f %.5f\n' %
                    (i_obj, r_icrs, d_icrs, np.degrees(ra), np.degrees(dec),
                    dmxx, dmyy, xxf, yyf, camxx, camyy))

    print('ra %.17f' % obs.pointingRA)
    print('dec %.17f' % obs.pointingDec)
    print('rotSky %.17f' % obs.rotSkyPos)
    print('mjd %.17f' % obs.mjd.TAI)
