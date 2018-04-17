import os
import numpy as np
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Bandpass, Sed
from lsst.sims.photUtils import PhysicalParameters
from lsst.sims.photUtils import PhotometricParameters

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    catsim_bp_dir = os.path.join(getPackageDir('throughputs'), 'imsim', 'goal')

    parser.add_argument('--phosim_dir', type=str,
                        help='Home directory for PhoSim',
                        default=os.path.join('/Users', 'danielsf', 'physics',
                                             'phosim_sfd'))

    parser.add_argument('--centroid_dir', type=str,
                        help='Directory under phosim_dir where centroid files are',
                        default = os.path.join('catsim_validation', 'hacked_filters'))

    parser.add_argument('--m1', type=str,
                        help='File to use for m1 reflectivity',
                        default=os.path.join(catsim_bp_dir, 'm1.dat'))

    parser.add_argument('--m2', type=str,
                        help='File to use for m2 reflectivity',
                        default=os.path.join(catsim_bp_dir, 'm2.dat'))

    parser.add_argument('--m3', type=str,
                        help='File to use for m3 reflectivity',
                        default=os.path.join(catsim_bp_dir, 'm3.dat'))

    parser.add_argument('--lens1', type=str,
                        help='File to use for lens1 transmission',
                        default=os.path.join(catsim_bp_dir, 'lens1.dat'))

    parser.add_argument('--lens2', type=str,
                        help='File to use for lens2 transmission',
                        default=os.path.join(catsim_bp_dir, 'lens2.dat'))

    parser.add_argument('--lens3', type=str,
                        help='File to use for lens3 transmission',
                        default=os.path.join(catsim_bp_dir, 'lens3.dat'))

    parser.add_argument('--det', type=str,
                        help='File to use for detector throughput (defaults to CatSim file)',
                        default=os.path.join(catsim_bp_dir, 'detector.dat'))

    parser.add_argument('--atm', type=str,
                        help='File to use for atmosphere throughput (defaults to airmass=1.2)',
                        default=os.path.join(catsim_bp_dir, 'atmos_std.dat'))

    parser.add_argument('--bp_dir', type=str,
                        help='Directory to search for filter_*.dat files',
                        default = catsim_bp_dir)

    args = parser.parse_args()
    phosim_dir = args.phosim_dir
    centroid_dir = os.path.join(phosim_dir, args.centroid_dir)
    m1_file = args.m1
    m2_file = args.m2
    m3_file = args.m3
    l1_file = args.lens1
    l2_file = args.lens2
    l3_file = args.lens3
    det_file = args.det
    atmos_file = args.atm
    bp_dir = args.bp_dir

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
    det = os.path.join(bp_dir, 'detector.dat')
    atmos = os.path.join(bp_dir, 'atmos_std.dat')

    componentList = [m1_file, m2_file, m3_file, l1_file, l2_file, l3_file]
    if det_file.lower() != 'none':
        print('adding detector')
        componentList.append(det_file)
    if atmos_file.lower() != 'none':
        print('adding astmosphere')
        componentList.append(atmos_file)

    optics_bp = Bandpass()
    optics_bp.readThroughputList(componentList=componentList)

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
