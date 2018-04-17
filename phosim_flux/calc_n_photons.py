import os
import numpy as np

_LSST_STACK_INSTALLED = True
try:
    from lsst.utils import getPackageDir
    from lsst.sims.photUtils import Bandpass, Sed
    from lsst.sims.photUtils import PhotometricParameters
except ImportError:
    _LSST_STACK_INSTALLED = False

import argparse

class PhysicalParameters(object):
    """
    A class to store physical constants and other immutable parameters
    used by the sims_photUtils code
    """

    def __init__(self):
        self.lightspeed = 299792458.0      # speed of light, = 2.9979e8 m/s
        self.planck = 6.626068e-27        # planck's constant, = 6.626068e-27 ergs*seconds
        self.nm2m = 1.00e-9               # nanometers to meters conversion = 1e-9 m/nm
        self.ergsetc2jansky = 1.00e23     # erg/cm2/s/Hz to Jansky units (fnu)


def get_sed_normalization(mag, sed_wav, sed_flambda):
    imsim_wav = np.arange(200.0, 700.0, 1.0)
    imsim_sb = np.zeros(len(imsim_wav))
    imsim_sb[300] = 1.0
    phi_integrand = imsim_sb/imsim_wav
    integral = 0.5*((phi_integrand[1:]+phi_integrand[:-1])
                    *(imsim_wav[1:]-imsim_wav[:-1])).sum()

    imsim_phi = phi_integrand/integral

    phys_params = PhysicalParameters()
    sed_fnu = sed_flambda*sed_wav*sed_wav
    sed_fnu *= phys_params.nm2m/phys_params.lightspeed
    sed_fnu *= phys_params.ergsetc2jansky

    sed_fnu = np.interp(imsim_wav, sed_wav, sed_fnu)

    current_flux = 0.5*((sed_fnu[1:]*imsim_phi[1:]+sed_fnu[:-1]*imsim_phi[:-1])
                        *(imsim_wav[1:]-imsim_wav[:-1])).sum()

    current_mag = -2.5*np.log10(current_flux) + 2.5*np.log10(3631)
    dmag = mag-current_mag
    fnorm = np.power(10.0, -0.4*dmag)
    return fnorm

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    try:
        catsim_bp_dir = os.path.join(getPackageDir('throughputs'), 'imsim', 'goal')
    except:
        catsim_bp_dir = '.'

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

    throughput_dtype = np.dtype([('wav_nm', float), ('throughput', float)])
    sed_dtype = np.dtype([('wav_nm', float), ('flambda', float)])

    phosim_ct_dtype = np.dtype([('id', int), ('phot', float), ('x', float), ('y', float)])

    sed_name = os.path.join(phosim_dir, 'data', 'SEDs',
                            'flatSED', 'sed_flat_short.txt.gz')

    np_sed = np.genfromtxt(sed_name, dtype=sed_dtype)
    np_fnorm = get_sed_normalization(21.0, np_sed['wav_nm'], np_sed['flambda'])
    np_sed['flambda'] *= np_fnorm

    if _LSST_STACK_INSTALLED:
        phot_params = PhotometricParameters(nexp=1, exptime=30.0)
        spec = Sed()
        spec.readSED_flambda(sed_name)
        imsim_bp = Bandpass()
        imsim_bp.imsimBandpass()
        fnorm = spec.calcFluxNorm(21.0, imsim_bp)
        spec.multiplyFluxNorm(fnorm)

    componentList = [m1_file, m2_file, m3_file, l1_file, l2_file, l3_file]
    if det_file.lower() != 'none':
        print('adding detector')
        componentList.append(det_file)
    if atmos_file.lower() != 'none':
        print('adding astmosphere')
        componentList.append(atmos_file)

    np_component_list = []
    for file_name in componentList:
        print('multiplying %s' % file_name)
        data = np.genfromtxt(file_name, dtype=throughput_dtype)
        np_component_list.append(data)

    if _LSST_STACK_INSTALLED:
        optics_bp = Bandpass()
        optics_bp.readThroughputList(componentList=componentList)

    for i_filter, bp_name in enumerate('ugrizy'):

        phosim_data = np.genfromtxt(os.path.join(centroid_dir,
                                                 'centroid_lsst_e_230_f%d_R22_S11_E000.txt' % i_filter),
                                    dtype=phosim_ct_dtype, skip_header=1)

        phosim_truth = np.median(phosim_data['phot'])

        np_filter = np.genfromtxt(os.path.join(bp_dir, 'filter_%s.dat' % bp_name),
                                  dtype=throughput_dtype)

        for np_component in np_component_list:
            interped_throughput = np.interp(np_filter['wav_nm'],
                                            np_component['wav_nm'],
                                            np_component['throughput'])

            np_filter['throughput'] *= interped_throughput


        if _LSST_STACK_INSTALLED:
            filter_bp = Bandpass()
            filter_bp.readThroughput(os.path.join(bp_dir, 'filter_%s.dat' % bp_name))

            wav, sb = optics_bp.multiplyThroughputs(filter_bp.wavelen, filter_bp.sb)
            bp = Bandpass(wavelen=wav, sb=sb)


        flambda = np.interp(np_filter['wav_nm'], np_sed['wav_nm'], np_sed['flambda'])
        phys_params = PhysicalParameters()
        phot = flambda*np_filter['wav_nm']/(phys_params.planck*phys_params.lightspeed*1.0e9)

        integral = 0.5*((phot[1:]*np_filter['throughput'][1:]+phot[:-1]*np_filter['throughput'][:-1])
                        *(np_filter['wav_nm'][1:]-np_filter['wav_nm'][:-1])).sum()

        effarea = np.pi*(6.423*100.0/2.0)**2
        exptime = 30.0
        nexp = 1.0

        integral *= effarea*exptime*nexp

        print('\n%s' % bp_name)
        print('catsim_counts (by hand) %e' % integral)

        if _LSST_STACK_INSTALLED:
            adu = spec.calcADU(bp, phot_params)
            print('catsim_counts (with sims_photUtils) %e' % (adu*phot_params.gain))
            if bp_name != 'y':
                assert np.abs(adu*phot_params.gain-integral) < 0.01*integral

        print('phosim_counts %e' % phosim_truth)
        print('catsim_counts/phosim_counts %e' % (integral/phosim_truth))
