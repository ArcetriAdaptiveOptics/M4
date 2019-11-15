'''
@author: cs
'''

from scipy import ndimage
import numpy as np
from m4.ground.configuration import Configuration
from m4.ground.zernikeGenerator import ZernikeGenerator
from m4.flattening import Flattenig

class ZernikeCommand():

    def __init__(self, segment_roi):
        self._roi = segment_roi
        self._pupilXYRadius = Configuration.M4_PUPIL_XYRADIUS
        self._zg = ZernikeGenerator(2*self._pupilXYRadius[2])
        self._segmentImaDiameter = 512

    def zernikeCommandTest(self, zernike_coeff_array, segment_number):
        surface_map = self._createZernikeModesOnM4(zernike_coeff_array)
        tetha, r = self._moveSegmentView(segment_number)

    def _createZernikeModesOnM4(self, zernike_coeff_array):
        ''' creao una superficie circolare con il diametro in pixel dello specchio'''
        zernike_surface_map = 0.0
        first_zern_mode_index = 2
        last_zern_mode_index = 2 + len(zernike_coeff_array)
        index_zernike_modes = np.arange(first_zern_mode_index, last_zern_mode_index)
        zd = self._zg.getZernikeDict(index_zernike_modes)

        for i in index_zernike_modes:
            zernike_surface_map = zernike_surface_map + \
                                    zernike_coeff_array[i-2] * zd[i]

        return zernike_surface_map


    def _cropImage(self, theta, r, zernike_surface_map):
        centrex = np.int(self._pupilXYRadius[0] + r * np.cos(theta))
        centrey = np.int(self._pupilXYRadius[1] + r * np.sin(theta))
        segment_radius = self._segmentImaDiameter // 2
        cropped_image = zernike_surface_map[centrey - segment_radius:
                                            centrey + segment_radius,
                                            centrex - segment_radius:
                                            centrex + segment_radius]
        rotate_image = ndimage.rotate(cropped_image, -theta - 90)
        final_image = np.ma.masked_array(rotate_image.data, mask= self._roi)

        return cropped_image, rotate_image, final_image


    def _moveSegmentView(self, segment_number):
        tetha = Configuration.REFERENCE_ANGLE * segment_number
        r = Configuration.SEGMENT_DISTANCE_FROM_CENTRE
        return tetha, r

        