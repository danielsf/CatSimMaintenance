import numpy as np
import os

__all__ = ["PhoSimPixelTransformer"]


def load_focal_plane():
    phosim_dir = os.path.join("/Users/danielsf/physics/phosim_release")
    phosim_file = os.path.join(phosim_dir, "data", "lsst", "focalplanelayout.txt")

    chip_data = {}
    with open(phosim_file, 'r') as input_file:
        for line in input_file:
            if line[0] == '#':
                continue
            params = line.strip().split()
            mangled_name = params[0]
            name = mangled_name[0]+':'+mangled_name[1]+','+mangled_name[2]+' '
            name += mangled_name[4]+':'+mangled_name[5]+','+mangled_name[6]
            local_dict = {}
            local_dict['x'] = float(params[1])
            local_dict['y'] = float(params[2])
            local_dict['p_size'] = float(params[3])
            local_dict['n_x'] = int(params[4])
            local_dict['n_y'] = int(params[5])
            local_dict['rot'] = float(params[10])
            assert np.abs(float(params[11]))<1.0e-10
            assert np.abs(float(params[12]))<1.0e-10
            local_dict['dx'] = float(params[13])
            local_dict['dy'] = float(params[14])
            assert np.abs(float(params[15]))<1.0e-10
            chip_data[name] = local_dict

    return chip_data


class PhoSimPixelTransformer(object):

    def __init__(self):
        self._chip_data = load_focal_plane()

    def _chip_center_mm(self, chipName):
        chip= self._chip_data[chipName]
        x0 = 0.001*chip['x']+chip['dx']
        y0 = 0.001*chip['y']+chip['dy']
        return x0, y0

    def _chip_rot_matrix(self, chipName):
        chip = self._chip_data[chipName]

        theta = np.radians(chip['rot'])
        print('theta %e' % chip['rot'])
        cc = np.cos(theta)
        ss = np.sin(theta)
        rotMat = np.array([[cc, -ss], [ss, cc]])
        return rotMat


    def mmFromPix(self, xpix, ypix, chipName):
        chip = self._chip_data[chipName]
        x0, y0 = self._chip_center_mm(chipName)
        dp = 0.001*chip['p_size']

        # rotate (1,0) and (0,1) so that x_vec and y_vec
        # reflect actual orientation of the pixel grid
        rotMat = self._chip_rot_matrix(chipName)

        x_pix_vec = np.array([1.0, 0.0])
        y_pix_vec = np.array([0.0, 1.0])
        x_pix_vec = np.dot(rotMat, x_pix_vec)
        y_pix_vec = np.dot(rotMat, y_pix_vec)

        mm_vec = np.array([x0, y0])
        mm_vec += (xpix-chip['n_x']/2)*dp*x_pix_vec
        mm_vec += (ypix-chip['n_y']/2)*dp*y_pix_vec

        return mm_vec[0], mm_vec[1]

    def pixFromMM(self, xmm, ymm, chipName):
        """
        xmm -- x coord in millimeters
        ymm -- y coord in millimeters
        chipName = like R:2,2 S:1,1
        """
        x0, y0 = self._chip_center_mm(chipName)

        chip = self._chip_data[chipName]
        dp = 0.001*chip['p_size']
        xpix0 = chip['n_x']/2
        ypix0 = chip['n_y']/2

        rotMat = self._chip_rot_matrix(chipName)

        x_pix_vec = np.array([1.0, 0.0])
        y_pix_vec = np.array([0.0, 1.0])
        x_pix_vec = np.dot(rotMat, x_pix_vec)
        y_pix_vec = np.dot(rotMat, y_pix_vec)

        d_mm_vec = np.array([xmm-x0, ymm-y0])

        dx = np.dot(x_pix_vec, d_mm_vec)
        dy = np.dot(y_pix_vec, d_mm_vec)

        xpix = xpix0 + dx/dp
        ypix = ypix0 + dy/dp

        return xpix, ypix
