"""
Like plot_displacement.py, but it uses obs_lsstSim to do all
coordinate transformations
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.GalSimInterface import LSSTCameraWrapper
from lsst.afw.cameraGeom import PIXELS, FOCAL_PLANE, FIELD_ANGLE
from lsst.afw.cameraGeom import SCIENCE
import lsst.afw.geom as afwGeom

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--obs', type=int, default=None)
parser.add_argument('--phosim_dir', type=str, default='phosim_output')
parser.add_argument('--out_dir', type=str, default='figs')
parser.add_argument('--catsim_dir', type=str, default='catalogs')

args = parser.parse_args()
if args.obs is None:
    raise RuntimeError("must specify obs")

phosim_dir = args.phosim_dir
out_dir = args.out_dir
assert os.path.isdir(phosim_dir)
phosim_file_list = os.listdir(phosim_dir)
for file_name in phosim_file_list:
    if 'centroid' in file_name:
        params = file_name.split('_')
        if int(params[3]) == args.obs:
            filter_name = params[4]
            break


camera = LsstSimMapper().camera
camera_wrapper = LSSTCameraWrapper()

det_name_list = []
for det in camera:
    if det.getType() != SCIENCE:
        continue
    det_name_list.append(det.getName())

phosim_dtype = np.dtype([('id', int), ('phot', int),
                         ('x', float), ('y', float)])

phosim_xmm = None
phosim_ymm = None
phosim_id = None
for det_name in det_name_list:
    detector = camera[det_name]
    m_name = det_name.replace(':','').replace(',','').replace(' ','_')
    file_name = 'centroid_lsst_e_%d_%s_%s_E000.txt' % (args.obs, filter_name, m_name)
    full_name = os.path.join(phosim_dir, file_name)

    if not os.path.exists(full_name):
        raise RuntimeError('%s does not exist' % full_name)

    local_data = np.genfromtxt(full_name, dtype=phosim_dtype, comments=None,
                               skip_header=1)

    x_dm, y_dm = camera_wrapper.dmPixFromCameraPix(local_data['x'], local_data['y'], det_name)

    pixels_to_focal = detector.getTransform(PIXELS, FOCAL_PLANE)

    x_mm = []
    y_mm = []
    for xx, yy in zip(x_dm, y_dm):
        focal = pixels_to_focal.applyForward(afwGeom.Point2D(xx, yy))
        x_mm.append(focal.getX())
        y_mm.append(focal.getY())
    
    x_mm = np.array(x_mm)
    y_mm = np.array(y_mm)

    if phosim_xmm is None:
        phosim_xmm = x_mm
        phosim_ymm = y_mm
        phosim_id = local_data['id']
    else:
        phosim_xmm = np.append(phosim_xmm, x_mm)
        phosim_ymm = np.append(phosim_ymm, y_mm)
        phosim_id = np.append(phosim_id, local_data['id'])


catsim_dtype = np.dtype([('id', int), ('xmm', float), ('ymm', float),
                         ('xpix', float), ('ypix', float)])

catsim_file = os.path.join(args.catsim_dir, 'star_predicted_%d.txt' % args.obs)
if not os.path.exists(catsim_file):
    raise RuntimeError('%s does not exist' % catsim_file)

catsim_data = np.genfromtxt(catsim_file, dtype=catsim_dtype)
catsim_id = catsim_data['id']
catsim_xmm = catsim_data['xmm']
catsim_ymm = catsim_data['ymm']

sorted_dex = np.argsort(catsim_id)
catsim_id = catsim_id[sorted_dex]
catsim_xmm = catsim_xmm[sorted_dex]
catsim_ymm = catsim_ymm[sorted_dex]

sorted_dex = np.argsort(phosim_id)
phosim_id = phosim_id[sorted_dex]
phosim_xmm = phosim_xmm[sorted_dex]
phosim_ymm = phosim_ymm[sorted_dex]

assert np.array_equal(phosim_id, catsim_id)

dx = phosim_xmm-catsim_xmm
dy = phosim_ymm-catsim_ymm

dd = np.sqrt(dx**2+dy**2)
print('dd %e %e %e' % (dd.min(),np.median(dd),dd.max()))

plt.figsize = (30,30)
plt.quiver(catsim_xmm, catsim_ymm, dx, dy)
plt.title('min %.2e med %.2e max %.2e' %
          (dd.min(), np.median(dd), dd.max()))
plt.savefig(os.path.join(out_dir,'catsim_to_phosim_%d.png') % args.obs)

dist_arr = np.sqrt(catsim_xmm**2+catsim_ymm**2)
disp_arr = np.sqrt(dx**2+dy**2)

with open(os.path.join(out_dir, 'catsim_to_phosim_%d_cat.txt' % args.obs), 'w') as out_file:
    out_file.write('# id catsim_xmm, catsim_ymm, phosim_xmm, phosim_ymm, dx, dy, total_displacement\n')
    for ii, x, y, px, py, d_x, d_y, ss in \
    zip(phosim_id, catsim_xmm, catsim_ymm, phosim_xmm, phosim_ymm, dx, dy, disp_arr):
        out_file.write('%d %.16e %.16e %.16e %.16e %e %e %e\n' % (ii, x, y, px, py, d_x, d_y, ss))
