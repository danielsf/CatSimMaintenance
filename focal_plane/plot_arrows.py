import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--chip', type=str, default=None)
args = parser.parse_args()
if args.chip is None:
    raise RuntimeError("must specify chip")

phosim_dir = os.path.join("/Users", "danielsf", "physics", "phosim_release")
assert os.path.isdir(phosim_dir)
centroid_file = os.path.join(phosim_dir, 'placement_grid',
                             'centroid_lsst_e_230_f2_%s_E000.txt' % args.chip)
assert os.path.exists(centroid_file)

phosim_dtype = np.dtype([('id', int), ('phot', int),
                         ('x', float), ('y', float)])

phosim_data = np.genfromtxt(centroid_file, dtype=phosim_dtype,
                            comments=None, skip_header=1)

catsim_dtype = np.dtype([('id', int), ('x', float), ('y', float),
                         ('xmm', float), ('ymm', float)])

catsim_file = os.path.join('raw_catalogs', 'star_cross_predicted_230_%s.txt' % args.chip)
assert os.path.exists(catsim_file)

catsim_data = np.genfromtxt(catsim_file, dtype=catsim_dtype)

print(len(phosim_data),len(catsim_data))

phosim_data = phosim_data[:-1]
assert len(catsim_data) == len(phosim_data)

assert np.array_equal(phosim_data['id'], catsim_data['id'])

plt.figsize = (30,30)

dx = phosim_data['x']-catsim_data['x']
dy = phosim_data['y']-catsim_data['y']

dd = np.sqrt(dx**2+dy**2)
print('dd %e %e %e' % (dd.min(),np.median(dd),dd.max()))

plt.quiver(catsim_data['x'], catsim_data['y'], dx, dy)
plt.savefig('catsim_to_phosim_%s.png' % args.chip)
