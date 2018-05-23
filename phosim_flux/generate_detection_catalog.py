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
                                          altitude=(50.0, 90.0),
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

xmm = rr*np.cos(theta)
ymm = rr*np.sin(theta)

xp, yp = pupilCoordsFromFocalPlaneCoordsLSST(xmm, ymm, band='r')
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
mag_norm = np.zeros(len(ra));
for ii in range(len(ra)):
    if ii%2==0:
        mag_norm[ii] = rng.random_sample()*2.0+17.0
    else:
        mag_norm[ii] =rng.random_sample()*4.0+18.5

out_dir = 'detection_catalogs_2'
sed_dir = 'cartoonSED'

sed_candidates = []
pwr = [1.0, -1.0, 2.0, -3.0, 3.0]
for ii in range(5):
    sed_name = os.path.join(sed_dir, 'sed_%d.txt' % ii)
    sed_candidates.append(sed_name)
    with open(sed_name, 'w') as out_file:
        wav = np.arange(200.0, 1400.0, 10.0)
        flambda = np.power(5.0*wav/wav.max(), pwr[ii])
        flambda = flambda/flambda.max()
        for ww, ff in zip(wav, flambda):
            out_file.write('%.5f %.5f\n' % (ww, ff))

        ss = Sed(wavelen=wav, flambda=flambda)
        mags = hw_dict.magListForSed(ss)
        for i_m in range(len(mags)-1):
            print('%.4e ' % (mags[i_m]-mags[i_m+1]))
        print('\n')

sed_names = rng.choice(sed_candidates, size=len(ra), replace=True)

opsim_meta_data = obs_root.OpsimMetaData

ref_made = False
ref_handle = None

mjd_grid = np.arange(59580.0, 59580.0+3652.5, 100.0)
mjd_good = []
for mjd in mjd_grid:
    obs_dummy = ObservationMetaData(mjd=mjd)
    alt, az, pa = altAzPaFromRaDec(obs_root.pointingRA, obs_root.pointingDec, obs_dummy)
    if alt>40.0:
        mjd_good.append(mjd)
    if len(mjd_good)>20:
        break


rot_sky_pos = rng.random_sample(len(mjd_good))*360.0
for i_mjd, mjd in enumerate(mjd_good):
    for i_filter in range(6):
        for all_same in (True, False):
            obs = ObservationMetaData(pointingRA=obs_root.pointingRA,
                                      pointingDec=obs_root.pointingDec,
                                      rotSkyPos = rot_sky_pos[i_mjd],
                                      mjd=mjd_good[i_mjd],
                                      site=obs_root.site,
                                      bandpassName='ugrizy'[i_filter])

            obs.OpsimMetaData = opsim_meta_data
            obs.OpsimMetaData['rawSeeing'] = rng.random_sample()*0.15+0.6

            if all_same:
                out_name = os.path.join(out_dir, 'instcat_%d_%d_same_sed.txt' % (i_mjd, i_filter))
            else:
                out_name = os.path.join(out_dir, 'instcat_%d_%d_diff_sed.txt' % (i_mjd, i_filter))

            ref_name = os.path.join(out_dir, 'phosim_ref_cat.txt')

            if not ref_made:
                ref_handle = open(ref_name, 'w')

            if all_same:
                obshistid = 1000
            else:
                obshistid = 2000
            obshistid += i_filter*100+i_mjd
            phosim_header['obshistid'] = obshistid

            if ref_handle is not None:
                ref_handle.write('# uniqueId, raJ2000, decJ2000, u, g, r, i, z, y, isresolved, isvariable\n')
            with open(out_name, 'w') as instcat_handle:
                write_phoSim_header(obs, instcat_handle, phosim_header)
                for ii, (rr, dd, mm) in enumerate(zip(ra, dec, mag_norm)):
                    instcat_handle.write('object %d ' % (ii+1))
                    instcat_handle.write('%.17f %.17f %.6f ' % (rr, dd, mm))
                    if ii%2==0 or all_same:
                        instcat_handle.write('flatSED/sed_flat_short.txt.gz ' )
                        if ii%2 == 0:
                            if ref_handle is not None:
                                ref_handle.write('%d, %.17f, %.17f, ' % (ii+1, rr, dd))
                                ref_handle.write('%.6f, %.6f, %.6f, %.6f, %.6f, %.6f, ' %
                                                 (mm, mm, mm, mm, mm, mm))
                                ref_handle.write('0, 0\n')
                    else:
                        instcat_handle.write('%s ' % sed_names[ii])
                    instcat_handle.write('0 0 0 0 0 0 point none none\n')

                if ref_handle is not None:
                    ref_handle.close()
                    ref_handle = None
