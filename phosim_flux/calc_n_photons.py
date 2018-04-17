import os
import numpy as np
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Bandpass, Sed
from lsst.sims.photUtils import PhysicalParameters
from lsst.sims.photUtils import PhotometricParameters

import argparse

if __name__ == "__main__":

    # directory where PhoSim lives
    phosim_dir = os.path.join('/Users', 'danielsf', 'physics',
                              'phosim_sfd')

    # directory where the centroid files are
    centroid_dir = os.path.join(phosim_dir, 'catsim_validation', 'hacked')

    phosim_ct_dtype = np.dtype([('id', int), ('phot', float), ('x', float), ('y', float)])

    phot_params = PhotometricParameters(nexp=1, exptime=30.0)

    spec = Sed()
    sed_name = os.path.join(phosim_dir, 'data', 'SEDs',
                            'flatSED', 'sed_flat_short.txt.gz')

    spec.readSED_flambda(sed_name)

    imsim_bp = Bandpass()
    imsim_bp.imsimBandpass()

    fnorm = spec.calcFluxNorm(21.0, imsim_bp)
    spec.multiplyFluxNorm(fnorm)

    bp_dir = os.path.join(getPackageDir('throughputs'), 'imsim', 'goal')
    m1 = os.path.join('hacked_throughputs', 'm1.dat')
    m2 = os.path.join('hacked_throughputs', 'm2.dat')
    m3 = os.path.join('hacked_throughputs', 'm3.dat')
    l1 = os.path.join('hacked_throughputs', 'lens1.dat')
    l2 = os.path.join('hacked_throughputs', 'lens2.dat')
    l3 = os.path.join('hacked_throughputs', 'lens3.dat')
    det = os.path.join(bp_dir, 'detector.dat')
    atmos = os.path.join(bp_dir, 'atmos_std.dat')

    optics_bp = Bandpass()
    optics_bp.readThroughputList(componentList=[m1, m2, m3, l1, l2, l3])

    phosim_bp_dtype = np.dtype([('angle', float), ('wav_micron', float), ('transmission', float),
                                ('reflection', float)])

    for i_filter, bp_name in enumerate('ugrizy'):

        phosim_data = np.genfromtxt(os.path.join(centroid_dir,
                                                 'centroid_lsst_e_230_f%d_R22_S11_E000.txt' % i_filter),
                                    dtype=phosim_ct_dtype, skip_header=1)

        phosim_truth = np.median(phosim_data['phot'])

        filter_bp = Bandpass()
        filter_bp.readThroughput(os.path.join(bp_dir, 'filter_%s.dat' % bp_name))
    

        wav, sb = optics_bp.multiplyThroughputs(filter_bp.wavelen, filter_bp.sb)
        bp = Bandpass(wavelen=wav, sb=sb)

        wav, flambda = spec.resampleSED(wavelen=spec.wavelen,
                                        flux=spec.flambda,
                                        wavelen_match=bp.wavelen)

        phys_params = PhysicalParameters()

        phot = flambda*wav/(phys_params.planck*phys_params.lightspeed*1.0e9)

        dlambda = wav[1]-wav[0]

        integral = 0.5*((phot[1:]*bp.sb[1:]+phot[:-1]*bp.sb[:-1])*(wav[1:]-wav[:-1])).sum()
        integral *= phot_params.effarea*phot_params.exptime*phot_params.nexp

        print('\n%s' % bp_name)
        print('catsim_counts %e' % integral)
        adu = spec.calcADU(bp, phot_params)
        assert np.abs(adu*phot_params.gain-integral) < 0.01*integral
        print('phosim_counts %e' % phosim_truth)
        print('catsim_counts/phosim_counts %e' % (integral/phosim_truth))
