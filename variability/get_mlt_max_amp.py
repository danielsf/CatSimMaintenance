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


def get_dust_grid():

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

    out_name = 'test_mlt_dmag.txt'
    if os.path.exists(out_name):
        raise RuntimeError("%s exists" % out_name)

    dust_grid = get_dust_grid()

    lc_data_dir = os.path.join(getPackageDir('sims_data'),
                               'catUtilsData')

    mdwarf_name = 'mdwarf_flare_light_curves_171012.npz'

    mag_name_list = ('u', 'g', 'r', 'i', 'z', 'y')
    flavor_to_int = {}
    flavor_to_int['late_active']  = 1
    flavor_to_int['late_inactive'] = 1
    flavor_to_int['mid_active'] = 2
    flavor_to_int['mid_inactive'] = 3
    flavor_to_int['early_active'] = 4
    flavor_to_int['early_inactive'] = 5

    lc_to_int = {}
    for i_lc in range(4):
        for flavor in ('late_active', 'late_inactive', 'early_active', 'early_inactive',
                       'mid_active', 'mid_inactive'):

            lc_to_int['%s_%d' % (flavor, i_lc)] = flavor_to_int[flavor]*10+i_lc

    lc_data = np.load(os.path.join(lc_data_dir, mdwarf_name))
    dflux_dict = {}
    for flavor in ('late_active', 'early_inactive', 'early_active',
                   'mid_active', 'mid_inactive'):
        for i_lc in range(4):
            lc_name = '%s_%d' % (flavor, i_lc)
            dflux_dict[lc_to_int[lc_name]] = {}
            for mag_name in mag_name_list:
                col_name = '%s_%d_%s' % (flavor, i_lc, mag_name)
                dflux = lc_data[col_name]
                dflux_dict[lc_to_int[lc_name]][mag_name] = np.abs(dflux).max()

    print('loaded dflux dict')
    ct_lt = 0
    ct_ge = 0
    cutoff=0.05
    table_tag = '1160'
    dtype=np.dtype([('htmid', int), ('simobjid', int),
                    ('umag', float), ('gmag', float), ('rmag', float),
                    ('imag', float), ('zmag', float), ('ymag', float),
                    ('sdssr', float),
                    ('varParamStr', str, 200), ('parallax', float),
                    ('ebv', float)])

    query = 'SELECT '
    query += 'htmid, simobjid, umag, gmag, rmag, imag, zmag, ymag, sdssr, '
    query += 'varParamStr, parallax, ebv '
    query += 'FROM stars_mlt_part_%s ' % table_tag
    query += 'TABLESAMPLE(1 percent)'

    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

    results_iter = db.get_arbitrary_chunk_iterator(query, dtype=dtype, chunk_size=100000)

    dummy_sed = Sed()
    _au_to_parsec = 1.0/206265.0
    _cm_per_parsec = 3.08576e18


    ct_bad = 0
    ct_total = 0
    for i_chunk, chunk in enumerate(results_iter):
        parallax = radiansFromArcsec(0.001*chunk['parallax'])
        dd = _au_to_parsec/parallax
        sphere_area = 4.0*np.pi*np.power(dd*_cm_per_parsec, 2)

        flux_factor = 1.0/sphere_area

        dust_factor = {}
        for mag_name in mag_name_list:
            dust_factor[mag_name] = np.interp(chunk['ebv'], dust_grid['ebv'], dust_grid[mag_name])

        lc_int_dex = np.ones(len(chunk), dtype=int)
        dmag_max = np.zeros(len(chunk), dtype=float)
        mag_min = np.zeros(len(chunk), dtype=float)
        for i_star, star in enumerate(chunk):
            var_dict = json.loads(star['varParamStr'])
            lc = var_dict['p']['lc']

            lc_int_dex[i_star] = lc_to_int[lc]

        for lc_id in np.unique(lc_int_dex):
            valid_stars = np.where(lc_int_dex==lc_id)
            if len(valid_stars[0])==0:
                continue
            local_flux_factor = flux_factor[valid_stars]
            local_dmag_max = np.zeros(len(valid_stars[0]), dtype=float)
            local_mag_min = 41.0*np.ones(len(valid_stars[0]), dtype=float)
            for mag_name in mag_name_list:
                mag_0 = chunk['%smag' % mag_name][valid_stars]
                local_dust_factor = dust_factor[mag_name][valid_stars]
                flux_0 = dummy_sed.fluxFromMag(mag_0)
                dflux = dflux_dict[lc_id][mag_name]*local_dust_factor*local_flux_factor
                total_flux = flux_0 + dflux
                mag = dummy_sed.magFromFlux(total_flux)
                dmag = 2.5*np.log10(1.0+dflux/flux_0)
                local_dmag_max = np.where(dmag>local_dmag_max, dmag, local_dmag_max)
                local_mag_min = np.where(mag<local_mag_min, mag, local_mag_min)

            dmag_max[valid_stars] = local_dmag_max
            mag_min[valid_stars] = local_mag_min

        ct_total += len(chunk)
        with open(out_name,'a') as out_file:
            for i_star in range(len(chunk)):
                if mag_min[i_star]>24.89 or dmag_max[i_star]<0.001:
                    ct_bad += 1
                out_file.write('%d %d %e %e\n' %
                (chunk['htmid'][i_star], chunk['simobjid'][i_star],
                 dmag_max[i_star],mag_min[i_star]))
        print('ct_bad %d of %d' % (ct_bad, ct_total))
