"""
This script will get all of the sources from a given OBAFGK table and
find a rough approximation of their maximum magnitude variability
"""
import os
import numpy as np
import json
import time
import gc
import multiprocessing as mproc
from lsst.sims.catalogs.db import DBObject
from lsst.sims.catUtils.mixins import ParametrizedLightCurveMixin
from lsst.sims.catUtils.mixins import create_variability_cache


def get_table_mins(table_tag, dmag_dict,_out_dir):
    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

    query = 'SELECT '
    query += 'htmid, simobjid, varParamStr '
    query += 'FROM stars_obafgk_part_%s' % table_tag

    dtype = np.dtype([('htmid', int), ('simobjid', int),
                      ('varParamStr', str, 200)])

    data_iter = get_arbitrary_chunk_iterator(query, dtype=dtype, chunk_size=10000)
    with open(os.path.join(_out_dir,'dmag_%s.txt' % table_tag), 'w') as out_file:
        for chunk in data_iter:
            for star in chunk:
                param_dict = json.loads(star['varParamStr'])
                lc_id = param_dict['p']['lc']
                out_file.write('%d %d %d\n' % (star['htmid'], star['simobjid'], dmag_dict[lc_id]))

if __name__ == "__main__":

    _out_dir = os.path.join('/astro/store/pogo4/danielsf/dmag_max_data')

    if not os.path.isdir(_out_dir):
        raise RuntimeError('%s is not a dir' % _out_dir)

    plc = ParametrizedLightCurveMixin()
    variability_cache = create_variability_cache()
    plc.load_parametrized_light_curves(variability_cache=variability_cache)
    time_arr = np.arange(0.0, 3653.0, 0.1)
    dmag_dict = {}

    t_start = time.time()
    for lc_id in variability_cache['_PARAMETRIZED_LC_MODELS']:
        q_flux, d_flux = plc._calc_dflux(lc_id, time_arr,
                                         variability_cache=variability_cache)
        d_flux = np.where(d_flux>-q_flux, d_flux, -0.9*q_flux)
        try:
            dmag = 2.5*np.log10(1.0+d_flux/q_flux)
            dmag = np.where(np.abs(dmag)>0.0, dmag, 1.0e-20)
            dmag = np.where(np.abs(dmag)<1.0e10, dmag, 1.0e10)
            bad = np.where(np.logical_or(np.isnan(dmag),np.isinf(dmag)))
            assert len(bad[0]) == 0
            dmag_dict[lc_id] = int(np.ceil(np.log10(np.abs(dmag).max())))
        except:
            print(lc_id)
            print(d_flux.min(),d_flux.max(),q_flux)
            print(dmag.min())
            print(dmag.max())
            is_inf = np.where(np.isinf(np.log10(np.abs(dmag).max())))
            print(is_inf)
            print(dmag[is_inf])
            raise
        if len(dmag_dict) % 1000 == 0:
            elapsed = (time.time()-t_start)/3600.0
            per = elapsed/len(dmag_dict)
            total = len(variability_cache['_PARAMETRIZED_LC_MODELS'])*per
            print('%d took %.2e hrs; per %.2e total %.2e' %
                  (len(dmag_dict),elapsed,per,total))

    del variability_cache
    gc.collect()

    with open(os.path.join(_out_dir,'dmag_dict.txt'),'w') as out_file:
        for lc_id in dmag_dict:
            out_file.write('%d %d\n' % (lc_id, dmag_dict[lc_id]))

    print('loaded dmag_dict')

    p_list = []
    for tag in ('0870', '1200', '1100', '1160', '1180', '1220', '1250', '1400'):
        p = mproc.Process(target=get_table_mins,
                          args=(tag, dmag_dict, _out_dir))

        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    print('all done')
