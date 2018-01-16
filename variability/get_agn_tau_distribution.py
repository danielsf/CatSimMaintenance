"""
Get the distribution of tau_agn in the AGN varParamStrs assigned on fatboy
"""

from lsst.sims.catalogs.db import DBObject
import json
import numpy as np

if __name__ == "__main__":

    db =DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                 port=1433, driver='mssql+pymssql')

    dtype = np.dtype([('varParamStr', str, 400), ('z', float)])

    query = 'SELECT varParamStr, redshift FROM galaxy WHERE varParamStr IS NOT NULL'

    data_iter = db.get_arbitrary_chunk_iterator(query, chunk_size=10000,
                                                dtype=dtype)

    with open('agn_tau_distribution.txt', 'w') as out_file:
        for chunk in data_iter:
            for line in chunk:
                params =  json.loads(line['varParamStr'])
                tau = params['pars']['agn_tau']
                out_file.write('%e %e\n' % (line['z'], tau))


