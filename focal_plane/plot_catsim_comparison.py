import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np
import argparse

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap
from lsst.sims.catUtils.exampleCatalogDefinitions import write_phoSim_header
from lsst.sims.utils import Site

from lsst.sims.utils import observedFromICRS
from lsst.sims.catUtils.mixins import PhoSimAstrometryBase
from lsst.sims.coordUtils import pixelCoordsFromRaDecLSST
from lsst.sims.coordUtils import lsst_camera

from PhoSimTransform import PhoSimPixelTransformer
from PhoSimTransform import _rawPupilCoordsFromObserved
from PhoSimTransform import _rawObservedFromPupilCoords
from lsst.sims.coordUtils import focalPlaneCoordsFromPupilCoords
from lsst.sims.coordUtils import pupilCoordsFromFocalPlaneCoords
from lsst.sims.utils import pupilCoordsFromRaDec

from lsst.sims.utils import altAzPaFromRaDec

if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument('--obs', type=int, default=230)
    parser.add_argument('--out_dir', type=str, default='catalogs')
    args = parser.parse_args()

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

    alt, az, pa = altAzPaFromRaDec(obs.pointingRA, obs.pointingDec,
                                   obs)

    assert alt>25.0

    transformer = PhoSimPixelTransformer(perturbed=True)
    phosim_mixin = PhoSimAstrometryBase()

    ra_icrs = np.arange(obs.pointingRA-0.15, obs.pointingRA+0.15, 0.02)
    dec_icrs = np.arange(obs.pointingDec-0.15, obs.pointingDec+0.15, 0.02)

    coord_grid = np.meshgrid(ra_icrs, dec_icrs)
    ra_icrs = coord_grid[0].flatten()
    dec_icrs = coord_grid[1].flatten()

    (xpup_catsim,
     ypup_catsim) = pupilCoordsFromRaDec(ra_icrs, dec_icrs,
                                         obs_metadata=obs,
                                         epoch=2000.0,
                                         includeRefraction=False)

    (x_focal_catsim,
     y_focal_catsim) = focalPlaneCoordsFromPupilCoords(xpup_catsim,
                                                       ypup_catsim,
                                                       camera=lsst_camera())

    ra_obs, dec_obs = observedFromICRS(ra_icrs, dec_icrs, obs_metadata=obs,
                                       epoch=2000.0, includeRefraction=False)

    ra_obs_rad = np.radians(ra_obs)
    dec_obs_rad = np.radians(dec_obs)

    (ra_deprecessed_rad,
     dec_deprecessed_rad) = phosim_mixin._dePrecess(ra_obs_rad, dec_obs_rad, obs)

    (xpup_deprecessed,
     ypup_deprecessed) = _rawPupilCoordsFromObserved(ra_deprecessed_rad,
                                                     dec_deprecessed_rad,
                                                     obs._pointingRA,
                                                     obs._pointingDec,
                                                     obs._rotSkyPos)

    (x_focal_deprecessed,
     y_focal_deprecessed) = focalPlaneCoordsFromPupilCoords(xpup_deprecessed,
                                                            ypup_deprecessed,
                                                            camera=lsst_camera())

    plt.figsize=(30,30)
    dx = x_focal_deprecessed-x_focal_catsim
    dy = y_focal_deprecessed-y_focal_catsim
    plt.quiver(x_focal_catsim, y_focal_catsim,
               dx, dy)
    plt.savefig('catsim_v_deprecessed_%d.eps' % args.obs)
