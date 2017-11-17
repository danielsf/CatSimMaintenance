"""
This script will get all of the sources from a given MLT table and
find a rough approximation of their maximum magnitude variability
"""

import os
import numpy as np
import json
from lsst.utils import getPackageDir
from lsst.sims.catalogs.db import DBObject
from lsst.sims.photUtils import Sed
from lsst.sims.utils import radiansFromArcsec


if __name__ == "__main__":

    lc_data_dir = os.path.join(getPackageDir('sims_data'),
                               'catUtilsData')

    mdwarf_name = 'mdwarf_flare_light_curves_171012.npz'

    mag_name_list = ('u', 'g', 'r', 'i', 'z', 'y')

    lc_data = np.load(os.path.join(lc_data_dir, mdwarf_name))
    dflux_dict = {}
    for flavor in ('late_active', 'early_inactive', 'early_active',
                   'mid_active', 'mid_inactive'):
        for i_lc in range(4):
            for mag_name in mag_name_list:
                col_name = '%s_%d_%s' % (flavor, i_lc, mag_name)
                dflux = lc_data[col_name]
                dflux_dict[col_name] = np.abs(dflux).max()

    print('loaded dflux dict')
    ct_lt = 0
    ct_ge = 0
    cutoff=0.05
    table_tag = '0870'
    dtype=np.dtype([('htmid', int), ('simobjid', int),
                    ('umag', float), ('gmag', float), ('rmag', float),
                    ('imag', float), ('zmag', float), ('ymag', float),
                    ('sdssr', float),
                    ('varParamStr', str, 200), ('parallax', float)])

    query = 'SELECT TOP 1000000'
    query += 'htmid, simobjid, umag, gmag, rmag, imag, zmag, ymag, sdssr, '
    query += 'varParamStr, parallax '
    query += 'FROM stars_mlt_part_%s' % table_tag

    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

    results_iter = db.get_arbitrary_chunk_iterator(query, dtype=dtype, chunk_size=10000)

    spec = Sed()
    _au_to_parsec = 1.0/206265.0
    _cm_per_parsec = 3.08576e18

    for chunk in results_iter:
        parallax = radiansFromArcsec(0.001*chunk['parallax'])
        dd = _au_to_parsec/parallax
        sphere_area = 4.0*np.pi*np.power(dd*_cm_per_parsec, 2)

        flux_factor = 1.0/sphere_area

        for i_star, star in enumerate(chunk):
            var_dict = json.loads(star['varParamStr'])
            lc = var_dict['p']['lc']
            if 'late' in lc:
                lc = lc.replace('inactive','active')
            lt = True

            for mag_name in mag_name_list:
                base_mag = star['%smag' % mag_name]
                base_flux = spec.fluxFromMag(base_mag)
                dmag = 2.5*np.log10(1.0 +dflux_dict['%s_%s' % (lc,mag_name)]*flux_factor[i_star]/base_flux)
                if np.abs(dmag) >= cutoff:
                    lt = False
                    break
            if lt:
                ct_lt += 1
            else:
                ct_ge += 1

            if dmag>9.0:
                print('lt %d ge %d %.2e %.2e %.2e %s' % (ct_lt, ct_ge, dmag, star['rmag'], star['sdssr'], lc))
