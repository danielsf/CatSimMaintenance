"""
Plot a quiver plot comparing pixel predictions from generate_grid_obs_lsstsim
with actual PhoSim outputs
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os

from lsst.sims.coordUtils import lsst_camera
from PhoSimTransform import load_focal_plane
from lsst.afw.cameraGeom import SCIENCE

if __name__ == "__main__":
    camera = lsst_camera()
    phosim_dir = 'phosim_full'
    assert os.path.isdir(phosim_dir)
    catsim_dir = 'full_catalogs'
    assert os.path.isdir(catsim_dir)

    phosim_dtype = np.dtype([('id', int), ('phot', float), ('x', float),
                             ('y', float)])

    phosim_id = None
    phosim_x = None
    phosim_y = None
    for det in camera:
        if det.getType()!=SCIENCE:
            continue
        chip_name = det.getName()
        chip_name_m = chip_name.replace(':','').replace(',','').replace(' ','_')
        file_name = 'centroid_lsst_e_230_f2_%s_E000.txt' % chip_name_m
        full_name = os.path.join(phosim_dir, file_name)
        try:
            assert os.path.exists(full_name)
        except AssertionError:
            raise RuntimeError('%s does not exist' % full_name)
        local_data = np.genfromtxt(full_name, dtype=phosim_dtype, skip_header=1)
        if phosim_id is None:
            phosim_id = local_data['id']
            phosim_x = local_data['x']
            phosim_y = local_data['y']
        else:
            phosim_id = np.append(phosim_id, local_data['id'])
            phosim_x = np.append(phosim_x, local_data['x'])
            phosim_y = np.append(phosim_y, local_data['y'])
    sorted_dex = np.argsort(phosim_id)

    phosim_id = phosim_id[sorted_dex]
    phosim_x = phosim_x[sorted_dex]
    phosim_y = phosim_y[sorted_dex]

    catsim_dtype = np.dtype([('id', int), ('ra', float), ('dec', float),
                             ('ra_deprecessed', float), ('dec_deprecessed', float),
                             ('x_dm', float), ('y_dm', float),
                             ('x_f', float), ('y_f', float),
                             ('x_cam', float), ('y_cam', float)])

    catsim_data = np.genfromtxt(os.path.join(catsim_dir, 'star_predicted_230.txt'),
                                dtype=catsim_dtype)

    print('len ',len(catsim_data))

    np.testing.assert_array_equal(phosim_id, catsim_data['id'])

    plt.figsize = (30, 30)
    dx = phosim_x-catsim_data['x_cam']
    dy = phosim_y-catsim_data['y_cam']
    dd = np.sqrt(dx**2+dy**2)
    print(dx)
    print(dy)
    print(dd)
    print(phosim_x)
    print(phosim_y)
    print(catsim_data['x_cam'])
    print(catsim_data['y_cam'])
    plt.quiver(catsim_data['x_f']/0.01, catsim_data['y_f']/0.01,
               dx, dy)
    plt.title('%.2e %.2e %.2e' % (dd.min(), np.median(dd), dd.max()))
    plt.savefig('full_focal_plan.eps')
