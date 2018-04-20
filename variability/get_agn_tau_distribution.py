"""
Get the distribution of AGN variability parameters assigned on fatboy
"""

from lsst.sims.catalogs.db import DBObject
import json
import numpy as np

if __name__ == "__main__":

    db =DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                 port=1433, driver='mssql+pymssql')

    dtype = np.dtype([('varParamStr', str, 400), ('z', float), ('magNorm', float)])

    query = 'SELECT varParamStr, redshift, magnorm_agn FROM galaxy WHERE varParamStr IS NOT NULL'

    data_iter = db.get_arbitrary_chunk_iterator(query, chunk_size=10000,
                                                dtype=dtype)

    with open('agn_variability_distribution.txt', 'w') as out_file:
        out_file.write('# z magnorm tau sfu sfg sfr sfi sfz sfy\n')
        for chunk in data_iter:
            for line in chunk:
                params =  json.loads(line['varParamStr'])
                tau = params['pars']['agn_tau']
                sfu = params['pars']['agn_sfu']
                sfg = params['pars']['agn_sfg']
                sfr = params['pars']['agn_sfr']
                sfi = params['pars']['agn_sfi']
                sfz = params['pars']['agn_sfz']
                sfy = params['pars']['agn_sfy']
                out_file.write('%e %e %e %e %e %e %e %e %e\n' %
                               (line['z'], line['magNorm'], tau,
                               sfu,sfg,sfr,sfi,sfz,sfy))


