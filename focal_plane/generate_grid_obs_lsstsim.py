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

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--obs', type=int, default=230)
    parser.add_argument('--chip', type=str, default='R:2,2 S:1,1')
    args = parser.parse_args()

    chip_m = args.chip.replace(':','').replace(',','').replace(' ','_')

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

    xpix_in = np.arange(200.0, 3800.0, 300.0)
    ypix_in = np.arange(200.0, 3800.0, 300.0)
    pix_grid = np.meshgrid(xpix_in, ypix_in)
    cam_xpix_in = pix_grid[0].flatten()
    cam_ypix_in = pix_grid[1].flatten()

    camera_wrapper = LSSTCameraWrapper()
    xpix_in, ypix_in = camera_wrapper.dmPixFromCameraPix(cam_xpix_in,
                                                         cam_ypix_in,
                                                         args.chip)

    ra_icrs, dec_icrs = raDecFromPixelCoords(xpix_in, ypix_in, args.chip,
                                             camera=lsst_camera(),
                                             obs_metadata=obs)

    ra_obs, dec_obs = observedFromICRS(ra_icrs, dec_icrs, obs_metadata=obs,
                                       epoch=2000.0, includeRefraction=False)

    de_precessor = PhoSimAstrometryBase()
    (ra_deprecessed,
     dec_deprecessed) = de_precessor._dePrecess(np.radians(ra_obs),
                                                np.radians(dec_obs),
                                                obs)

    phosim_header_map = copy.deepcopy(DefaultPhoSimHeaderMap)
    phosim_header_map['nsnap'] = 1
    phosim_header_map['exptime'] = 30.0

    with open('full_catalogs/test_cat_%s_%d.txt' % (chip_m, args.obs), 'w') as out_file:
        write_phoSim_header(obs, out_file, phosim_header_map)
        for i_obj, (ra, dec) in enumerate(zip(ra_deprecessed, dec_deprecessed)):
            out_file.write('object %d ' % (i_obj+1))
            out_file.write('%.17f %.17f ' % (np.degrees(ra), np.degrees(dec)))
            out_file.write('21.0 flatSED/sed_flat_short.txt.gz 0 0 0 0 0 0 point none CCM 0.03380581 3.1\n')

    with open('full_catalogs/input_pix_%s_%d.txt' % (chip_m,args.obs), 'w') as out_file:
        out_file.write('# i_obj x_cam y_cam\n')
        for i_obj, (xx, yy) in enumerate(zip(cam_xpix_in, cam_ypix_in)):
            out_file.write('%d %e %e\n' % (i_obj+1, xx, yy))
