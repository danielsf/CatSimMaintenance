import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

from PhoSimTransform import PhoSimPixelTransformer
from lsst.sims.coordUtils import lsst_camera
from lsst.afw.cameraGeom import SCIENCE

import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--obs', type=int, default=None)
parser.add_argument('--out_dir', type=str, default='figs')
parser.add_argument('--phosim_dir', type=str, default='phosim_output')
args = parser.parse_args()
if args.obs is None:
    raise RuntimeError("must specify obs")

phosim_dtype = np.dtype([('id', int), ('phot', int),
                         ('x', float), ('y', float)])



camera = lsst_camera()
det_name_list = []
for det in camera:
    if det.getType() != SCIENCE:
        continue
    det_name_list.append(det.getName())

det_name_list.sort()
det_name_map = {}
for i_det, det_name in enumerate(det_name_list):
    det_name_map[det_name] = i_det

coordinate_converter = PhoSimPixelTransformer()

position_dict = {}
for i_filter in range(6):
    chip_grid = []
    chip_dex_grid = []
    phosim_x = None
    phosim_y = None
    phosim_id = None
    for det_name in det_name_list:
        m_name = det_name.replace(':','').replace(',','').replace(' ','_')
        file_name = 'centroid_lsst_e_%d_f%d_%s_E000.txt' % (args.obs, i_filter, m_name)
        full_name = os.path.join(args.phosim_dir, file_name)

        local_data = np.genfromtxt(full_name, dtype=phosim_dtype, comments=None,
                                   skip_header=1)

        if phosim_x is None:
            phosim_x = local_data['x']
            phosim_y = local_data['y']
            phosim_id = local_data['id']
        else:
            phosim_x = np.append(phosim_x, local_data['x'])
            phosim_y = np.append(phosim_y, local_data['y'])
            phosim_id = np.append(phosim_id, local_data['id'])

        chip_grid += [det_name]*len(local_data)
        chip_dex_grid += [det_name_map[det_name]]*len(local_data)

    phosim_xmm = []
    phosim_ymm = []
    for name, xx, yy in zip(chip_grid, phosim_x, phosim_y):
        xmm, ymm = coordinate_converter.mmFromPix(xx, yy, name)
        phosim_xmm.append(xmm)
        phosim_ymm.append(ymm)

    phosim_xmm = np.array(phosim_xmm)
    phosim_ymm = np.array(phosim_ymm)
    chip_dex_grid = np.array(chip_dex_grid)

    position_dict[i_filter] = {}
    position_dict[i_filter]['id'] = phosim_id
    position_dict[i_filter]['xmm'] = phosim_xmm
    position_dict[i_filter]['ymm'] = phosim_ymm
    position_dict[i_filter]['xpix'] = phosim_x
    position_dict[i_filter]['ypix'] = phosim_y
    position_dict[i_filter]['chip'] = chip_dex_grid

for i_filter in range(6):
    assert np.array_equal(position_dict[i_filter]['id'], position_dict[0]['id'])

with open(os.path.join(args.out_dir, 'position_by_filter.txt'), 'w') as out_file:
    out_file.write('# xmm_u ymm_u xmm_g ymm_g xmm_r ymm_r...\n')
    for ii in range(len(position_dict[0]['id'])):
        for i_filter in range(6):
            out_file.write('%.6f %.6f ' % (position_dict[i_filter]['xmm'][ii],
                                           position_dict[i_filter]['ymm'][ii]))
        out_file.write('\n')

for i_filter in range(6):
    if i_filter==2:
        continue

    filt_name= ('u', 'g', 'r', 'i', 'z', 'y')[i_filter]
    plt.figsize=(300, 300)
    dx = position_dict[i_filter]['xmm'] - position_dict[2]['xmm']
    dy = position_dict[i_filter]['ymm'] - position_dict[2]['ymm']
    plt.quiver(position_dict[2]['xmm'], position_dict[2]['ymm'],
               dx, dy)

    plt.xlabel('xmm (r filter)')
    plt.ylabel('ymm (r filter)')

    plt.savefig(os.path.join(args.out_dir, 'r_to_%s_mm.eps' % filt_name))
    plt.close()


for i_filter in range(6):
    if i_filter==2:
        continue
    filt_name= ('u', 'g', 'r', 'i', 'z', 'y')[i_filter]

    for chip_name in ('R:0,3 S:1,1', 'R:4,3 S:1,1',
                      'R:0,1 S:1,1', 'R:4,1 S:1,1'):

        mangled_name = chip_name.replace(':','').replace(',','').replace(' ','_')

        chip_dex = det_name_map[chip_name]
        valid_obj = np.where(position_dict[2]['chip'] == chip_dex)
        x_r = position_dict[2]['xpix'][valid_obj]
        y_r = position_dict[2]['ypix'][valid_obj]
        id_r = position_dict[2]['id'][valid_obj]
        x_f = position_dict[i_filter]['xpix'][valid_obj]
        y_f = position_dict[i_filter]['ypix'][valid_obj]
        id_f = position_dict[i_filter]['id'][valid_obj]
        assert np.array_equal(id_r, id_f)

        plt.figsize=(300, 300)
        dx = x_f-x_r
        dy = y_f-y_r
        plt.quiver(x_r, y_r, dx, dy)
        plt.xlabel('x pixels')
        plt.ylabel('y pixels')
        plt.savefig(os.path.join(args.out_dir, 'r_to_%s_%s_pix.eps' % (filt_name,mangled_name)))
        plt.close()
