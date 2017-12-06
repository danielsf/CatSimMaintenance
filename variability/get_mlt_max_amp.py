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
from lsst.sims.photUtils import BandpassDict


def get_mlt_grid():

    # Construct a look-up table to determine the factor
    # by which to multiply the flares' flux to account for
    # dust as a function of E(B-V).  Recall that we are
    # modeling all MLT flares as 9000K blackbodies.

    lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

    ebv_grid = np.arange(0.0, 7.01, 0.01)
    bb_wavelen = np.arange(200.0, 1500.0, 0.1)
    hc_over_k = 1.4387e7  # nm*K
    temp = 9000.0  # black body temperature in Kelvin
    exp_arg = hc_over_k/(temp*bb_wavelen)
    exp_term = 1.0/(np.exp(exp_arg) - 1.0)
    ln_exp_term = np.log(exp_term)

    # Blackbody f_lambda function;
    # discard normalizing factors; we only care about finding the
    # ratio of fluxes between the case with dust extinction and
    # the case without dust extinction
    log_bb_flambda = -5.0*np.log(bb_wavelen) + ln_exp_term
    bb_flambda = np.exp(log_bb_flambda)
    bb_sed = Sed(wavelen=bb_wavelen, flambda=bb_flambda)

    base_fluxes = lsstBandpassDict.fluxListForSed(bb_sed)

    a_x, b_x = bb_sed.setupCCMab()
    _mlt_dust_lookup = {}
    _mlt_dust_lookup['ebv'] = ebv_grid
    list_of_bp = lsstBandpassDict.keys()
    for bp in list_of_bp:
        _mlt_dust_lookup[bp] = np.zeros(len(ebv_grid))
    for iebv, ebv_val in enumerate(ebv_grid):
        wv, fl = bb_sed.addCCMDust(a_x, b_x,
                                   ebv=ebv_val,
                                   wavelen=bb_wavelen,
                                   flambda=bb_flambda)

        dusty_bb = Sed(wavelen=wv, flambda=fl)
        dusty_fluxes = lsstBandpassDict.fluxListForSed(dusty_bb)
        for ibp, bp in enumerate(list_of_bp):
            _mlt_dust_lookup[bp][iebv] = dusty_fluxes[ibp]/base_fluxes[ibp]

    return _mlt_dust_lookup



if __name__ == "__main__":

    mlt_grid = get_mlt_grid()

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
                    ('varParamStr', str, 200), ('parallax', float),
                    ('ebv', float)])

    query = 'SELECT TOP 1000000'
    query += 'htmid, simobjid, umag, gmag, rmag, imag, zmag, ymag, sdssr, '
    query += 'varParamStr, parallax, ebv '
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
