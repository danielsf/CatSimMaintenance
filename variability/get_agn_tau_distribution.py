"""
Get the distribution of AGN variability parameters assigned on fatboy
"""

from lsst.sims.catalogs.db import DBObject
import json
import numpy as np
import os

from lsst.utils import getPackageDir
from lsst.sims.photUtils import BandpassDict, Sed, Bandpass

if __name__ == "__main__":

    db =DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                 port=1433, driver='mssql+pymssql')

    dtype = np.dtype([('varParamStr', str, 400), ('z', float), ('magNorm', float)])

    query = 'SELECT varParamStr, redshift, magnorm_agn FROM galaxy WHERE varParamStr IS NOT NULL'

    data_iter = db.get_arbitrary_chunk_iterator(query, chunk_size=10000,
                                                dtype=dtype)



    sed_file = os.path.join(getPackageDir('sims_sed_library'),
                            'agnSED', 'agn.spec.gz')

    raw_agn_sed = Sed()
    raw_agn_sed.readSED_flambda(sed_file)

    imsim_bp = Bandpass()
    imsim_bp.imsimBandpass()

    bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

    with open('agn_variability_distribution.txt', 'w') as out_file:
        out_file.write('# z u g r i z y tau sfu sfg sfr sfi sfz sfy\n')
        for chunk in data_iter:
            for line in chunk:
                agn_sed = Sed(wavelen=raw_agn_sed.wavelen,
                              flambda=raw_agn_sed.flambda)

                fnorm = agn_sed.calcFluxNorm(line['magNorm'], imsim_bp)
                agn_sed.multiplyFluxNorm(fnorm)
                agn_sed.redshiftSED(line['z'], dimming=True)
                mags = bp_dict.magListForSed(agn_sed)
                params =  json.loads(line['varParamStr'])
                tau = params['pars']['agn_tau']
                sfu = params['pars']['agn_sfu']
                sfg = params['pars']['agn_sfg']
                sfr = params['pars']['agn_sfr']
                sfi = params['pars']['agn_sfi']
                sfz = params['pars']['agn_sfz']
                sfy = params['pars']['agn_sfy']
                out_file.write('%e %e %e %e %e %e %e %e %e %e %e %e %e %e\n' %
                               (line['z'],
                                mags[0], mags[1], mags[2], mags[3],
                                mags[4], mags[5], tau,
                                sfu,sfg,sfr,sfi,sfz,sfy))


