'''
@author: cs
'''

import numpy as np
from m4.ground.configuration import Configuration
from m4.ground.zernikeGenerator import ZernikeGenerator
from m4.ground.zernikeMask import CircularMask


class ZernikeOnM4():
    """
    Class for the generation of Zernike modes in relation to the deformable mirror.

    HOW TO USE IT:
    from m4.utils.zernike_on_m_4 import ZernikeOnM4
    zOnM4= ZernikeOnM4()
    """

    def __init__(self):
        """The constructor """
        self._pupilXYRadius = Configuration.PARABOLA_PUPIL_XYRADIUS
        self._zg = ZernikeGenerator(2*self._pupilXYRadius[2])

    def getPupilCenterAndRadiusInIFCoords(self):
        return self._pupilXYRadius

    def setPupilCenterAndRadiusInIFCoords(self, centerX, centerY, radius):
        self._pupilXYRadius = np.array([centerX, centerY, radius])
        self._zg = ZernikeGenerator(2*radius)


    def zernikeFit(self, img, zernike_mode):
        '''
        arg:
            img = numpy masked array
            zernike_mode = vector of Zernike modes to remove

        return:
            a = vector containing Zernike modes amplitudes
            mat = interaction matrix for the masked area
        '''
        mat = np.zeros((img.compressed().shape[0], zernike_mode.size))
        for i in range(0, zernike_mode.size):
            z = self._zg.getZernike(zernike_mode[i])
            aa = np.ma.masked_array(z, mask=img.mask)
            mat.T[i] = aa.compressed()

        self._mat = mat
        inv = np.linalg.pinv(mat)
        a = np.dot(inv, img.compressed())
        return a, mat

    def zernikeSurface(self, surface_zernike_coeff_array, ima_mask,
                       mat, index=None):
        '''
        args:
            surface_zernike_coeff_array = vector containing the amplitudes
                                        of the Zernike modes
            ima_mask = the mask area in which we want to rebuilding
                    the surfaces
            mat = interaction matrix for the masked area
            index = vector containing the index number of interaction matrix
                    that we want to use
        returns:
            surf = reconstructed surface
        '''
        zernike_surface_map = None

        if index is None:
            zernike_surface_map = np.dot(mat, surface_zernike_coeff_array)
        else:
            for i in range(len(index)):
                k = index[i]
                zernike_surface = np.dot(mat[:, k],
                                         surface_zernike_coeff_array[i])
                if zernike_surface_map is None:
                    zernike_surface_map = zernike_surface
                else:
                    zernike_surface_map = zernike_surface_map + zernike_surface

        mask = np.invert(ima_mask)
        surf = np.ma.masked_array(np.zeros((2*self._pupilXYRadius[2],
                                            2*self._pupilXYRadius[2])),
                                  mask=mask)
        surf[mask] = zernike_surface_map
        return surf


    def zernikeToDMCommand(self, surface_map, an):
        '''
        Args:
            surface_map =
            an = analyzer delle funzioni d'influenza zonali
        Returns:
            zernike_cmd = Command (numpy.array)
        '''
        self._an = an
        cmask = self._an.getMasterMask()
        zernike_mask_on_if = CircularMask(self._an.getIFShape(),
                                          self._pupilXYRadius[2],
                                          [self._pupilXYRadius[1],
                                           self._pupilXYRadius[0]]).zernikeMask()
        self._an.setDetectorMask(zernike_mask_on_if | cmask)
        rec = self._an.getReconstructor()

        zernike_cmd = np.dot(rec, surface_map.compressed())

        return zernike_cmd
