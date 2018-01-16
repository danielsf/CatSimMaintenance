import numpy as np
from lsst.sims.catUtils.baseCatalogModels import StarObj
from lsst.sims.utils import getAllTrixels, levelFromHtmid
from lsst.sims.utils import halfSpaceFromRaDec, trixelFromHtmid, findHtmid
from lsst.sims.catUtils.utils import StellarAlertDBObjMixin

class StarObjHTMID(StarObj, StellarAlertDBObjMixin):
    pass


def get_kurucz_phys(sed_name):
    """
    Read in the name of a kurucz SED file.  Return it's
    T_eff, metallicity, log(surface gravity)
    """
    if sed_name[1] == 'p':
        metallicity_sgn = 1.0
    elif sed_name[1] == 'm':
        metallicity_sgn = -1.0
    else:
        raise RuntimeError('Cannot parse metallicity sign of %s' % sed_name)

    new_name = sed_name.replace('.','_').split('_')

    metallicity = 0.1*metallicity_sgn*float(new_name[0][2:])

    teff = float(new_name[-1])

    logg = 0.1*np.float(new_name[3][1:])

    return teff, metallicity, logg


def get_wd_phys(sed_name):
    """
    Read in the name of a white dwarf SED,
    return its T_eff, metallicity (which we don't actually have),
    and log(surface gravity)
    """
    new_name = sed_name.replace('.','_').split('_')
    teff = float(new_name[-1])
    if new_name[1]!='He':
        logg = 0.1*float(new_name[2])
    else:
        logg = 0.1*float(new_name[3])

    return teff, -999.0, logg


def get_mlt_phys(sed_name):
    """
    Read in the name of an M/L/T dwarf SED and return
    its T_eff, metallicity, and log(surface gravity)
    """

    new_name = sed_name.replace('+','-').replace('a','-').split('-')

    logg_sgn_dex = len(new_name[0])

    if sed_name[logg_sgn_dex] == '-':
        logg_sgn = 1.0
    elif sed_name[logg_sgn_dex] == '+':
        logg_sgn = -1.0
    else:
        raise RuntimeError('Cannot get logg_sgn for %s' % sed_name)

    metallicity_sgn_dex = len(new_name[0]) + len(new_name[1]) + 1

    if sed_name[metallicity_sgn_dex] == '-':
        metallicity_sgn = -1.0
    elif sed_name[metallicity_sgn_dex] == '+':
        metallicity_sgn = 1.0
    else:
        raise RuntimeError('Cannot get metallicity_sgn for %s' % sed_name)

    teff = 100.0*float(new_name[0][3:])
    metallicity = metallicity_sgn*float(new_name[2])
    logg = logg_sgn*float(new_name[1])

    return teff, metallicity, logg


def get_physical_characteristics(sed_name):
    """
    Read in the name of an SED file.
    Return (in this order) Teff, metallicity (FeH), log(g)
    """
    sed_name = sed_name.strip()
    if sed_name.endswith('.gz'):
        sed_name = sed_name.replace('.gz','')
    if sed_name.endswith('.txt.'):
        sed_name = sed_name.replace('.txt','')

    if not hasattr(get_physical_characteristics, 'teff_dict'):
        get_physical_characteristics.teff_dict = {}
        get_physical_characteristics.logg_dict = {}
        get_physical_characteristics.metal_dict = {}

    if sed_name in get_physical_characteristics.teff_dict:
        return (get_physical_characteristics.teff_dict[sed_name],
                get_physical_characteristics.metal_dict[sed_name],
                get_physical_characteristics.logg_dict[sed_name])

    if sed_name.startswith('bergeron'):
        sub_dir = 'wDs'
    elif sed_name.startswith('lte'):
        sub_dir = 'mlt'
    elif sed_name[0] == 'k':
        sub_dir = 'kurucz'
    else:
        raise RuntimeError("Do not understand name %s" % sed_name)



    if 'kurucz' in sub_dir:
        tt,mm,gg =  get_kurucz_phys(sed_name)
    elif sub_dir == 'wDs':
        tt,mm,gg = get_wd_phys(sed_name)
    elif sub_dir == 'mlt':
        tt,mm,gg = get_mlt_phys(sed_name)
    else:
        raise RuntimeError('Do not know how to get '
                           'physical characteristics for '
                           'sub_dir %s' % sub_dir)

    get_physical_characteristics.teff_dict[sed_name] = tt
    get_physical_characteristics.metal_dict[sed_name] = mm
    get_physical_characteristics.logg_dict[sed_name] = gg

    return tt, mm, gg


if __name__ == "__main__":

    trixel_dict = getAllTrixels(6)

    hs1 = halfSpaceFromRaDec(0.0, 90.0, 91.3)
    hs2 = halfSpaceFromRaDec(0.0, -90.0, 91.3)
    for ra in range(0,360,20):
        htmid = findHtmid(ra,0.0,max_level=6)
        tx = trixelFromHtmid(htmid)
        assert hs1.contains_trixel(tx) != 'outside'
        assert hs2.contains_trixel(tx) != 'outside'


    valid_htmid = []
    n_6 = 0
    for htmid in trixel_dict:
        if levelFromHtmid(htmid) != 6:
            continue
        n_6 += 1
        tx = trixel_dict[htmid]
        if hs1.contains_trixel(tx) != 'outside':
            if hs2.contains_trixel(tx) != 'outside':
                valid_htmid.append(htmid)
    print(n_6,len(valid_htmid))

    db = StarObjHTMID(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                      port=1433, driver='mssql+pymssql')

    colnames = ['ra', 'decl', 'gal_l', 'gal_b',
                'mura', 'mudecl', 'parallax',
                'vrad', 'sedfilename',
                'umag', 'gmag', 'rmag', 'imag', 'zmag', 'ymag']

    constraint = 'rmag<=23.0 AND rmag>=13.0 AND '
    constraint += '(ra>=300.0 OR ra<=60.0) AND '
    constraint += 'decl>=-1.27 AND decl<=1.27'
    ct = 0
    with open('stripe_82_catalog.txt', 'w') as out_file:
        out_file.write('# Columns are:\n')
        out_file.write('# RA (degrees)\n')
        out_file.write('# Dec (degrees)\n')
        out_file.write('# galactic longitude (degrees)\n')
        out_file.write('# galactic lagitutde (degrees)\n')
        out_file.write('# proper motion RA (mas/yr)\n')
        out_file.write('# proper motion Dec (mas/yr)\n')
        out_file.write('# parallax (mas)\n')
        out_file.write('# radial velocity (km/s)\n')
        out_file.write('# Teff (Kelvin)\n')
        out_file.write('# FeH\n')
        out_file.write('# log(g)\n')
        out_file.write('# umag\n')
        out_file.write('# gmag\n')
        out_file.write('# rmag\n')
        out_file.write('# imag\n')
        out_file.write('# zmag\n')
        out_file.write('# ymag\n')

        for i_htmid, htmid in enumerate(valid_htmid):

            data_iter = db.query_columns_htmid(colnames=colnames, constraint=constraint,
                                               chunk_size=100000, htmid=htmid)


            for chunk in data_iter:
                ct += len(chunk)
                for star in chunk:
                    params = get_physical_characteristics(star['sedfilename'])
                    out_file.write('%.6f %.6f ' % (star['ra'], star['decl']))
                    out_file.write('%.6f %.6f ' % (star['gal_l'], star['gal_b']))
                    out_file.write('%e %e %e ' % (star['mura'], star['mudecl'], star['parallax']))
                    out_file.write('%e %e %e %e ' % (star['vrad'], params[0], params[1], params[2]))
                    out_file.write('%.4f %.4f %.4f ' % (star['umag'], star['gmag'], star['rmag']))
                    out_file.write('%.4f %.4f %.4f ' % (star['imag'], star['zmag'], star['ymag']))
                    out_file.write('\n')
                print('%d stars -- i_htmid %d' % (ct, i_htmid))
