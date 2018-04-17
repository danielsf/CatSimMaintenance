"""
Like generate_grid_catalog.py, but using obs_lsstSim to map from
pixels to pupil and focal plane coordinates
"""

import os
import numpy as np
import copy

from PhoSimTransform import _rawObservedFromPupilCoords
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.utils import Site

import lsst.afw.geom as afwGeom

from lsst.afw.cameraGeom import SCIENCE
from lsst.afw.cameraGeom import FIELD_ANGLE, FOCAL_PLANE, PIXELS

from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap
from lsst.sims.catUtils.exampleCatalogDefinitions import write_phoSim_header

import argparse

if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument('--obs', type=int, default=230)
    parser.add_argument('--out_dir', type=str, default='catalogs')
    parser.add_argument('--chip', type=str, default=None)
    parser.add_argument('--random', type=str, default='False')

    args = parser.parse_args()
    if args.random.lower()[0] == 't':
        rng = np.random.RandomState(9912312)
    else:
        rng = None

    opsimdb = os.path.join('/Users', 'danielsf', 'physics', 'lsst_150412',
                           'Development', 'garage', 'OpSimData',
                           'minion_1016_sqlite.db')

    assert os.path.exists(opsimdb)
    obs_gen = ObservationMetaDataGenerator(database=opsimdb)
    obs_list = obs_gen.getObservationMetaData(obsHistID=args.obs)
    obs = obs_list[0]
    filter_name= obs.bandpass
    site_no_atm = Site(name="LSST",
                       pressure=0.0,
                       humidity=0.0)
    obs.site=site_no_atm
    assert np.abs(obs.site.pressure)<1.0e-6
    assert np.abs(obs.site.humidity)<1.0e-6
    print('pointing ',obs.pointingRA,obs.pointingDec)

    camera = LsstSimMapper().camera

    x_pix_arr = np.arange(700.0, 3301.0, 700.0)
    y_pix_arr = np.arange(700.0, 3301.0, 700.0)
    pix_grid = np.meshgrid(x_pix_arr, y_pix_arr)
    x_pix_arr = pix_grid[0].flatten()
    y_pix_arr = pix_grid[1].flatten()

    x_mm = []
    y_mm = []
    x_pup = []
    y_pup = []
    det_name = []
    id_grid = []

    det_name_list = []
    for det in camera:
        if det.getType() != SCIENCE:
            continue
        if args.chip is not None:
            if det.getName() != args.chip:
                continue
        det_name_list.append(det.getName())

    det_name_list.sort()

    focal_to_field = camera.getTransformMap().getTransform(FOCAL_PLANE,
                                                           FIELD_ANGLE)

    obj_id = 1
    for dn in det_name_list:
        detector = camera[dn]
        pixels_to_focal = detector.getTransform(PIXELS, FOCAL_PLANE)

        if rng is not None:
            x_pix_arr = rng.random_sample(25)*3500.0 + 300.0
            y_pix_arr = rng.random_sample(25)*3500.0 + 300.0

        for ii in range(len(x_pix_arr)):
            focal = pixels_to_focal.applyForward(afwGeom.Point2D(x_pix_arr[ii],
                                                                 y_pix_arr[ii]))
            field = focal_to_field.applyForward(focal)
            x_mm.append(focal.getX())
            y_mm.append(focal.getY())
            x_pup.append(field.getX())
            y_pup.append(field.getY())
            det_name.append(dn)
            id_grid.append(obj_id)
            obj_id += 1

    x_pup = np.array(x_pup)
    y_pup = np.array(y_pup)

    ra, dec = _rawObservedFromPupilCoords(x_pup, y_pup,
                                          obs._pointingRA, obs._pointingDec,
                                          obs._rotSkyPos)

    ra_grid = np.degrees(ra)
    dec_grid = np.degrees(dec)

    name_to_num ={'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5}

    phosim_header_map = copy.deepcopy(DefaultPhoSimHeaderMap)
    phosim_header_map['rawSeeing'] = ('rawSeeing', None)
    phosim_header_map['FWHMeff'] = ('FWHMeff', None)
    phosim_header_map['FWHMgeom'] = ('FWHMgeom',None)
    phosim_header_map['nsnap'] = 1
    phosim_header_map['vistime'] =30.0
    obs.OpsimMetaData['obsHistID'] = name_to_num[filter_name]

    with open(os.path.join(args.out_dir,'star_grid_%d.txt' % (name_to_num[filter_name])), 'w') as out_file:
        write_phoSim_header(obs, out_file, phosim_header_map)
        for (i_obj, ra, dec) in zip(id_grid, ra_grid, dec_grid):
            out_file.write('object %d ' % (i_obj))
            out_file.write('%.17f %.17f ' % (ra, dec))
            out_file.write('21.0 flatSED/sed_flat_short.txt.gz 0 0 0 0 0 0 point none CCM 0.03380581 3.1\n')

    with open(os.path.join(args.out_dir, 'star_predicted_%d.txt' % (name_to_num[filter_name])), 'w') as out_file:
        out_file.write('# ra %.17e dec %.17e rotSkyPos %.17e MJD(TAI) %.17e\n' %
                       (obs.pointingRA, obs.pointingDec, obs.rotSkyPos, obs.mjd.TAI))
        out_file.write('# id xmm ymm xpup ypup ra_obs dec_obs\n')
        for ii in range(len(x_mm)):
            out_file.write('%d %.17e %.17e %.17e %.17e %.17e %.17e\n' %
                           (id_grid[ii],
                            x_mm[ii], y_mm[ii],
                            x_pup[ii], y_pup[ii],
                            ra_grid[ii], dec_grid[ii]))

