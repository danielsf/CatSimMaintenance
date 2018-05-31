import os
import numpy as np
from lsst.utils import getPackageDir
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap
from lsst.sims.catUtils.exampleCatalogDefinitions import write_phoSim_header
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.coordUtils import pupilCoordsFromFocalPlaneCoordsLSST
from lsst.sims.utils import raDecFromPupilCoords
from lsst.sims.utils import angularSeparation
from lsst.sims.utils import Site
from lsst.sims.photUtils import BandpassDict, Sed
from lsst.sims.utils import altAzPaFromRaDec
from lsst.sims.utils import ObservationMetaData

bp_dict, hw_dict = BandpassDict.loadBandpassesFromFiles()

opsimdb = os.path.join('/Users/danielsf/physics/lsst_150412',
                       'Development', 'garage', 'OpSimData',
                       'minion_1016_sqlite.db')

obs_gen = ObservationMetaDataGenerator(opsimdb)
obs_list = obs_gen.getObservationMetaData(moonAlt=(-90.0, -50.0),
                                          altitude=(55.0, 57.0),
                                          fieldDec=(-10.0, 10.0))

assert len(obs_list) > 0
obs_root = obs_list[0]

obs_root.site = Site(name='LSST', pressure=0.0, humidity=0.0)

phosim_header = DefaultPhoSimHeaderMap
phosim_header['nsnap'] = 1
phosim_header['vistime'] = 30.0

galaxy_dir = os.path.join(getPackageDir('sims_sed_library'), 'galaxySED')
galaxy_sed_list = os.listdir(galaxy_dir)


magnorm_interp = {}

from lsst.sims.photUtils import PhotometricParameters, Bandpass
phot_params = PhotometricParameters(nexp=1, exptime=30.0, gain=1)
imsim = Bandpass()
imsim.imsimBandpass()

target=20000.0
for name in galaxy_sed_list:
    magnorm_interp[name] = {}
    magnorm_interp[name]['z'] = np.arange(0.0, 2.4, 0.2)
    magnorm_interp[name]['mag'] = np.zeros(len(magnorm_interp[name]['z']), dtype=float)
    for i_zz, zz in enumerate(magnorm_interp[name]['z']):
        ss = Sed()
        ss.readSED_flambda(os.path.join(galaxy_dir, name))
        mag = ss.calcMag(imsim)
        ss.redshiftSED(zz, dimming=True)
        cts = ss.calcADU(bp_dict['r'], photParams=phot_params)
        magnorm = mag -2.5*np.log10(target/cts)

        ss = Sed()
        ss.readSED_flambda(os.path.join(galaxy_dir, name))
        fnorm = ss.calcFluxNorm(magnorm, imsim)
        ss.multiplyFluxNorm(fnorm)
        ss.redshiftSED(zz, dimming=True)
        new_cts = ss.calcADU(bp_dict['r'], photParams=phot_params)
        d_ct = np.abs(new_cts-target)
        
        magnorm_interp[name]['mag'][i_zz] = magnorm
        
        if d_ct>500.0:
            raise RuntimeWarning('d_ct %e' % d_ct)

rng = np.random.RandomState(18113)
rr_chip = np.sqrt(2)*127.0*0.5/3.0 # radius of a chip in mm
n_sources = 200

rr = rng.random_sample(n_sources)*rr_chip
theta = rng.random_sample(n_sources)*2.0*np.pi

from lsst.sims.coordUtils import pupilCoordsFromPixelCoordsLSST
from lsst.sims.coordUtils import lsst_camera
from lsst.afw.cameraGeom import SCIENCE

xpix = np.arange(200.0, 3801.0, 200.0)
ypix = np.arange(200.0, 3801.0, 200.0)
mesh = np.meshgrid(xpix, ypix)
xpix_mesh =mesh[0].flatten()
ypix_mesh = mesh[1].flatten()

xpix = []
ypix = []
chip = []
for det in lsst_camera():
    if det.getType()!=SCIENCE:
        continue
    for xx, yy in zip(xpix_mesh, ypix_mesh):
        xpix.append(xx)
        ypix.append(yy)
        chip.append(det.getName())

xpix = np.array(xpix)
ypix = np.array(ypix)
chip = np.array(chip)
xp, yp = pupilCoordsFromPixelCoordsLSST(xpix, ypix, chipName=chip, band='r')
ra, dec = raDecFromPupilCoords(xp, yp, obs_metadata=obs_root,
                               includeRefraction=False)

out_dir = 'ratio_test/catalogs'

sed_candidates = os.listdir(galaxy_dir)

sed_names = rng.choice(sed_candidates, size=len(ra), replace=True)
redshift_vals = rng.random_sample(len(sed_names))*2.0

for i_filter in range(6):
    obs = ObservationMetaData(pointingRA=obs_root.pointingRA,
                              pointingDec=obs_root.pointingDec,
                              rotSkyPos=obs_root.rotSkyPos,
                              mjd=obs_root.mjd,
                              site=obs_root.site,
                              bandpassName='ugrizy'[i_filter])
    obs.OpsimMetaData = obs_root.OpsimMetaData
    obs.OpsimMetaData['sunalt'] = -90.0
    phosim_header['sunalt'] = -90.0
    phosim_header['camconfig'] = 1

    out_name = os.path.join(out_dir, 'instcat_ratio_grid_%d.txt' % (i_filter))
    with open(out_name, 'w') as instcat_handle:
        phosim_header['obshistid'] = 100+i_filter

        write_phoSim_header(obs, instcat_handle, phosim_header)
        obj_ct = 0
        redshift = 0.0
        for ii, (rr, dd, zz, ss) in enumerate(zip(ra, dec, redshift_vals, sed_names)):
            magnorm = np.interp(zz, magnorm_interp[ss]['z'], magnorm_interp[ss]['mag'])
            unqid = ii+1
            instcat_handle.write('object %d ' % unqid)
            instcat_handle.write('%.17f %.17f %.6f ' % (rr, dd, magnorm))
            instcat_handle.write('galaxySED/%s ' % ss)
            instcat_handle.write('%.6f 0 0 0 0 0 point none none\n' % zz)
