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

bp_dict, hw_dict = BandpassDict.loadBandpassesFromFiles()

opsimdb = os.path.join('/Users/danielsf/physics/lsst_150412',
                       'Development', 'garage', 'OpSimData',
                       'minion_1016_sqlite.db')

obs_gen = ObservationMetaDataGenerator(opsimdb)
obs_list = obs_gen.getObservationMetaData(moonAlt=(-90.0, -50.0),
                                          altitude=(50.0, 90.0),
                                          fieldDec=(-10.0, 10.0))

assert len(obs_list) > 0
obs = obs_list[0]

obs.site = Site(name='LSST', pressure=0.0, humidity=0.0)

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
ra, dec = raDecFromPupilCoords(xp, yp, obs_metadata=obs,
                               includeRefraction=False)

to_drop = []

fwhm = obs.OpsimMetaData['FWHMeff']
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
print(obs.OpsimMetaData['FWHMeff'])

to_drop.sort(reverse=True)
ra = list(ra)
dec = list(dec)
for ii in to_drop:
    ra.pop(ii)
    dec.pop(ii)

ra = np.array(ra)
dec = np.array(dec)
mag_norm = rng.random_sample(len(ra))*6.0+18.0

out_dir = 'detection_catalogs'
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

for i_filter in [0]:
    out_name = os.path.join(out_dir, 'phosim_flux_cat_%d.txt' % i_filter)
    phosim_header['obshistid'] = i_filter
    with open(out_name, 'w') as file_handle:
        write_phoSim_header(obs, file_handle, phosim_header)
        for ii, (rr, dd, mm) in enumerate(zip(ra, dec, mag_norm)):
            file_handle.write('object %d ' % (ii+1))
            file_handle.write('%.17f %.17f %.6f ' % (rr, dd, mm))
            if ii%2 == 0:
                file_handle.write('flatSED/sed_flat_short.txt.gz ' )
            else:
                file_handle.write('%s ' % sed_names[ii])
            file_handle.write('0 0 0 0 0 0 point none none\n')
