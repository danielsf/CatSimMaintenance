import numpy as np
from PhoSimTransform import PhoSimPixelTransformer
from lsst.afw.cameraGeom import SCIENCE
from lsst.sims.coordUtils import lsst_camera


if __name__ == "__main__":

    coord_converter = PhoSimPixelTransformer()
    camera = lsst_camera()
    rng = np.random.RandomState(213)
    det_name_list = []
    for det in camera:
        if det.getType() != SCIENCE:
            continue
        det_name_list.append(det.getName())

    n_test = 20
    test_name_list = rng.choice(det_name_list, size=n_test, replace=True)
    xpix_list = rng.random_sample(n_test)*4000.0
    ypix_list = rng.random_sample(n_test)*4000.0
    for xpix, ypix, test_name in zip(xpix_list, ypix_list, test_name_list):

        xmm, ymm = coord_converter.mmFromPix(xpix, ypix, test_name)
        xpix1, ypix1 = coord_converter.pixFromMM(xmm, ymm, test_name)
        dd = np.sqrt((xpix-xpix1)**2+(ypix-ypix1)**2)
        assert dd<1.0e-10
