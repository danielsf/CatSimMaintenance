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

control_sed = 'flatSED/sed_flat_short.txt.gz'


rng = np.random.RandomState(18113)
rr_chip = np.sqrt(2)*127.0*0.5/3.0 # radius of a chip in mm
n_sources = 200

rr = rng.random_sample(n_sources)*rr_chip
theta = rng.random_sample(n_sources)*2.0*np.pi

from lsst.sims.coordUtils import pupilCoordsFromPixelCoordsLSST
from lsst.sims.coordUtils import lsst_camera
from lsst.afw.cameraGeom import SCIENCE

xx = [1000.0, 1500.0, 2000.0, 2500.0, 3000.0]
yy = [2000.0, 2000.0, 2000.0, 2000.0, 2000.0]

xpix = []
ypix = []
chip = []
for det in lsst_camera():
    if det.getType() != SCIENCE:
        continue

    for x, y in zip(xx, yy):
        xpix.append(x)
        ypix.append(y)
        chip.append(det.getName())
        

xpix = np.array(xpix)
ypix = np.array(ypix)
chip = np.array(chip)
xp, yp = pupilCoordsFromPixelCoordsLSST(xpix, ypix, chipName=chip, band='r')
ra, dec = raDecFromPupilCoords(xp, yp, obs_metadata=obs_root,
                               includeRefraction=False)

to_drop = []

fwhm = obs_root.OpsimMetaData['FWHMeff']
min_dist = 10.0*fwhm/3600.0

for ii, (rr, dd) in enumerate(zip(ra, dec)):
    dd = angularSeparation(rr, dd, ra, dec)
    too_close = np.where(dd<min_dist)
    for i_too in too_close[0]:
        if i_too!=ii and i_too not in to_drop:
            to_drop.append(ii)
            break

print(len(ra))
print(len(to_drop))
print(obs_root.OpsimMetaData['FWHMeff'])

to_drop.sort(reverse=True)
ra = list(ra)
dec = list(dec)
for ii in to_drop:
    ra.pop(ii)
    dec.pop(ii)

ra = np.array(ra)
dec = np.array(dec)

magnorm = 21.0
redshift = 0.8

out_dir = 'galaxy_grid'

sed_dir = getPackageDir('sims_sed_library')
galaxy_dir = os.path.join(sed_dir, 'galaxySED')
sed_candidates = os.listdir(galaxy_dir)

print(len(sed_candidates))

sed_names = rng.choice(sed_candidates, size=len(ra), replace=False)

opsim_meta_data = obs_root.OpsimMetaData

ref_made = False
ref_handle = None

from lsst.sims.photUtils import Bandpass, PhotometricParameters
phot_params = PhotometricParameters(nexp=1, exptime=30.0)

imsim = Bandpass()
imsim.imsimBandpass()

for redshift in np.arange(0.0, 0.8, 0.2):
    for all_same in (True, False):

        ref_name = os.path.join(out_dir, 'galaxy_grid_ref_cat_%.1f.txt' % redshift)


        with open(ref_name,'w') as ref_handle:
            ref_handle.write('# uniqueId, raJ2000, decJ2000, u, g, r, i, z, y, isresolved, isvariable\n')
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

                out_name = os.path.join(out_dir, 'instcat_galaxy_grid_%.1f_%d.txt' % (redshift,i_filter))
                with open(out_name, 'w') as instcat_handle:
                    phosim_header['obshistid'] = i_filter*100+int(redshift/0.1)+1

                    write_phoSim_header(obs, instcat_handle, phosim_header)
                    for ii, (rr, dd) in enumerate(zip(ra, dec)):
                        unqid = ii+1
                        instcat_handle.write('object %d ' % unqid)
                        instcat_handle.write('%.17f %.17f %.6f ' % (rr, dd, magnorm))
                        instcat_handle.write('galaxySED/%s ' % sed_names[ii])
                        instcat_handle.write('%1.f 0 0 0 0 0 point none none\n' % redshift)
                        if i_filter==0:
                            ss = Sed()
                            ss.readSED_flambda(os.path.join(galaxy_dir, sed_names[ii]))
                            fnorm = ss.calcFluxNorm(magnorm, imsim)
                            ss.multiplyFluxNorm(fnorm)
                            ss.redshiftSED(redshift, dimming=True)
                            ref_handle.write('%d %.17f %.17f ' % (unqid, rr, dd))
                            for bp in bp_dict:
                                adu = ss.calcADU(bp_dict[bp], photParams=phot_params)
                                ref_handle.write('%e ' % adu)
                            ref_handle.write('\n')
