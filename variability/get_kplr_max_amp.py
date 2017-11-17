"""
This script will get all of the sources from a given OBAFGK table and
find a rough approximation of their maximum magnitude variability
"""
import os
import numpy as np
import json
from lsst.utils import getPackageDir
from lsst.sims.catalogs.db import DBObject
from lsst.sims.photUtils import Sed
from lsst.sims.utils import radiansFromArcsec
from lsst.sims.catUtils.mixins import ParametrizedLightCurveMixin

if __name__ == "__main__":

    plc = ParametrizedLightCurveMixin()
    plc.load_parametrized_light_curves()
    time_arr = np.arange(0.0, 3653.0, 0.1)
    dmag_dict = {}

    table_tag = '0870'
    dtype=np.dtype([('htmid', int), ('simobjid', int),
                    ('umag', float), ('gmag', float), ('rmag', float),
                    ('imag', float), ('zmag', float), ('ymag', float),
                    ('sdssr', float),
                    ('varParamStr', str, 200), ('parallax', float)])

    query = 'SELECT TOP 1000000'
    query += 'htmid, simobjid, umag, gmag, rmag, imag, zmag, ymag, sdssr, '
    query += 'varParamStr, parallax '
    query += 'FROM stars_obafgk_part_%s' % table_tag

    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

    results_iter = db.get_arbitrary_chunk_iterator(query, dtype=dtype, chunk_size=10000)

    cutoff= 0.05
    ct_lt = 0
    ct_ge = 0
    for chunk in results_iter:
        for star in chunk:
            var_dict = json.loads(star['varParamStr'])
            lc_id = var_dict['p']['lc']
            if lc_id not in dmag_dict:
                q_flux, d_flux = plc._calc_dflux(lc_id, time_arr)
                dmag = 2.5*np.log10(1.0+d_flux/q_flux)
                dmag_max = np.abs(dmag).max()
                dmag_dict[lc_id] = dmag_max
            dmag = dmag_dict[lc_id]
            if dmag<cutoff:
                ct_lt +=1
            else:
                ct_ge += 1
            print('lt %d ge %d -- %d' % (ct_lt,ct_ge,len(dmag_dict)))
