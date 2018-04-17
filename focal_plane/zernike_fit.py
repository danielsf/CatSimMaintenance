import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np
from lsst.obs.lsstSim import LsstSimMapper
from lsst.afw.cameraGeom import PIXELS, FOCAL_PLANE, SCIENCE
import lsst.afw.geom as afwGeom
from lsst.sims.GalSimInterface import LSSTCameraWrapper
from lsst.sims.utils import ZernikePolynomialGenerator

from lsst.sims.coordUtils import LsstZernikeFitter

import time

z_fitter = LsstZernikeFitter()
camera = LsstSimMapper().camera
camera_wrapper = LSSTCameraWrapper()

t_start = time.time()

z_gen = ZernikePolynomialGenerator()
n_grid = []
m_grid = []
for n in range(4):
    for m in range(-n, n+1, 2):
        n_grid.append(n)
        m_grid.append(m)

fig_dir = os.path.join('learning', 'figs')
catalog_dir = os.path.join('learning', 'catalogs')
phosim_dir = os.path.join('learning', 'output')
prediction_cat = os.path.join(catalog_dir, 'star_predicted_2.txt')
dtype = np.dtype([('id', int), ('xmm', float), ('ymm', float),
                  ('xpup', float), ('ypup', float),
                  ('raObs', float), ('decObs', float)])

catsim_data = np.genfromtxt(prediction_cat, dtype=dtype)

sorted_dex = np.argsort(catsim_data['id'])
catsim_data=catsim_data[sorted_dex]
for ii in range(len(catsim_data)):
    assert catsim_data['id'][ii] == (ii+1)

dtype = np.dtype([('id', int), ('phot', float),
                  ('xpix', float), ('ypix', float)])

catsim_radius = np.sqrt(catsim_data['xmm']**2 + catsim_data['ymm']**2)
rr = catsim_radius.max()
catsim_x = catsim_data['xmm']/rr
catsim_y = catsim_data['ymm']/rr

polynomials ={}
for n, m in zip(n_grid, m_grid):
    values = z_gen.evaluate_xy(catsim_x, catsim_y, n, m)
    polynomials[(n,m)] = values

poly_keys = list(polynomials.keys())

for i_filter in range(6):
    dx = np.zeros(len(catsim_data['xpup']), dtype=float)
    dy = np.zeros(len(catsim_data['ypup']), dtype=float)
    phosim_xmm = np.zeros(len(catsim_data['ypup']), dtype=float)
    phosim_ymm = np.zeros(len(catsim_data['ypup']), dtype=float)

    for det in camera:
        if det.getType() != SCIENCE:
            continue
        pixels_to_focal = det.getTransform(PIXELS, FOCAL_PLANE)
        det_name = det.getName()
        det_name_m = det_name.replace(':','').replace(',','').replace(' ','_')
        centroid_name = 'centroid_lsst_e_2_f%d_%s_E000.txt' % (i_filter, det_name_m)
        full_name = os.path.join(phosim_dir, centroid_name)
        phosim_data = np.genfromtxt(full_name, dtype=dtype, skip_header=1)
        xpix, ypix = camera_wrapper.dmPixFromCameraPix(phosim_data['xpix'], phosim_data['ypix'],
                                                       det_name)
        xmm = np.zeros(len(xpix), dtype=float)
        ymm = np.zeros(len(ypix), dtype=float)
        for ii in range(len(xpix)):
            focal_pt = pixels_to_focal.applyForward(afwGeom.Point2D(xpix[ii], ypix[ii]))
            xmm[ii] = focal_pt.getX()
            ymm[ii] = focal_pt.getY()
        dx[phosim_data['id']-1] = xmm - catsim_data['xmm'][phosim_data['id']-1]
        dy[phosim_data['id']-1] = ymm - catsim_data['ymm'][phosim_data['id']-1]
        phosim_xmm[phosim_data['id']-1] = xmm
        phosim_ymm[phosim_data['id']-1] = ymm

    b = np.array([(dx*polynomials[k]).sum() for k in poly_keys])
    m = np.array([[(polynomials[k1]*polynomials[k2]).sum() for k1 in poly_keys]
                  for k2 in poly_keys])

    alpha_x = np.linalg.solve(m, b)
    dx_fit = np.zeros(len(dx))
    for ii in range(len(alpha_x)):
        dx_fit += alpha_x[ii]*polynomials[poly_keys[ii]]

    #print('\n x')
    ddx = np.abs(dx-dx_fit)
    ddx_sorted = np.sort(ddx)
    n = len(ddx_sorted)
    #print(ddx_sorted[0], ddx_sorted[n//4], ddx_sorted[n//2], ddx_sorted[3*n//4])
    abs_dx = np.abs(dx)
    abs_dx_sorted = np.sort(abs_dx)
    #print(abs_dx_sorted[0], abs_dx_sorted[n//4], abs_dx_sorted[n//2], abs_dx_sorted[3*n//4])

    b = np.array([(dy*polynomials[k]).sum() for k in poly_keys])
    m = np.array([[(polynomials[k1]*polynomials[k2]).sum() for k1 in poly_keys]
                  for k2 in poly_keys])

    alpha_y = np.linalg.solve(m, b)
    dy_fit = np.zeros(len(dy))
    for ii in range(len(alpha_y)):
        dy_fit += alpha_y[ii]*polynomials[poly_keys[ii]]

    #print('\ny')

    ddy = np.abs(dy-dy_fit)
    ddy_sorted = np.sort(ddy)
    n = len(ddy_sorted)
    #print(ddy_sorted[0], ddy_sorted[n//4], ddy_sorted[n//2], ddy_sorted[3*n//4])
    abs_dy = np.abs(dy)
    abs_dy_sorted = np.sort(abs_dy)
    #print(abs_dy_sorted[0], abs_dy_sorted[n//4], abs_dy_sorted[n//2], abs_dy_sorted[3*n//4])

    disp = np.sqrt(dx**2+dy**2)
    new_x = catsim_data['xmm'] + dx_fit
    new_y = catsim_data['ymm'] + dy_fit
    new_disp = np.sqrt((phosim_xmm-new_x)**2+(phosim_ymm-new_y)**2)
    
    disp_sorted = np.sort(disp)/0.01
    new_disp_sorted = np.sort(new_disp)/0.01
    n = len(disp)
    print('\ntotal %d' % i_filter)
    print('new: %e %e %e %e' % (new_disp_sorted[n//4],
                                new_disp_sorted[n//2], new_disp_sorted[3*n//4],
                                new_disp_sorted[-1]))


    print('old: %e %e %e %e' % (disp_sorted[n//4],
                                disp_sorted[n//2], disp_sorted[3*n//4],
                                disp_sorted[-1]))

    dx_catsim, dy_catsim = z_fitter.dxdy(catsim_data['xmm'],
                                         catsim_data['ymm'],
                                         i_filter)

    np.testing.assert_array_almost_equal(dx_catsim, dx_fit, decimal=10)
    np.testing.assert_array_almost_equal(dy_catsim, dy_fit, decimal=10)

print('non camera took %.2e ' % (time.time()-t_start))
