'''
Autors
  - C. Selmi: written in 2019
'''

import numpy as np
from m4.configuration.ott_parameters import *
from m4.ground.zernikeGenerator import ZernikeGenerator
from m4.ground.zernikeMask import CircularMask
from m4.ground import geo
from m4.ground import zernike


class ZernikeOnM4():
    """
    Class for the generation of Zernike modes in relation to the deformable mirror.

    HOW TO USE IT::

        from m4.utils.zernike_on_m_4 import ZernikeOnM4
        zOnM4= ZernikeOnM4()
    """

    def __init__(self):
        """The constructor """
        self._pupilXYRadius = OttParameters.PARABOLA_PUPIL_XYRADIUS
        self._zg = ZernikeGenerator(2*self._pupilXYRadius[2])

    def getPupilCenterAndRadiusInIFCoords(self):
        '''
        Returns
        -------
        pupilXYRadius: numpy array
            vector containing x,y coordinates of center and radius
            used to create the Zernike
        '''
        return self._pupilXYRadius

    def setPupilCenterAndRadiusInIFCoords(self, centerX, centerY, radius):
        '''
        Parameters
        ----------
        centerX: x coordinate of the center
        centerY: y coordinate of the center
        radius: radius of the pupil
        '''
        self._pupilXYRadius = np.array([centerX, centerY, radius])
        self._zg = ZernikeGenerator(2*radius)

    def zernikeFit(self, img, zernike_index_vector):
        img1 = img.data
        mask = np.invert(img.mask).astype(int)
        x, y, r, xx, yy = geo.qpupil(mask)
        mm = (mask==1)
        coeff = zernike.surf_fit(xx[mm], yy[mm], img1[mm], zernike_index_vector)
        mat = zernike.getZernike(xx[mm], yy[mm], zernike_index_vector)
        return coeff, mat

    def zernikeSurface(self, img, zernike_index_vector):
        img1 = img.data
        mask = np.invert(img.mask).astype(int)
        x, y, r, xx, yy = geo.qpupil(mask)
        mm = (mask==1)
        aa = zernike.getZernike(xx[mm], yy[mm], zernike_index_vector)
        coeff = zernike.surf_fit(xx[mm], yy[mm], img1[mm], zernike_index_vector)
        zernike_surface = np.zeros((img.shape[0], img.shape[1]))
        zernike_surface[mm] = np.dot(aa, coeff)
        surf = np.ma.masked_array(zernike_surface, mask=img.mask)
        return surf





 ### Funzioni vecchie #####       

    def _zernikeFit(self, img, zernike_mode):
        '''
        Parameters
        ----------
            img: numpy masked array
            zernike_mode: numpy array
                    vector of Zernike modes to remove

        Returns
        -------
            a: numpy array
                vector containing Zernike modes amplitudes
            mat: numpy array
                interaction matrix for the masked area
        '''
        mat = np.zeros((img.compressed().shape[0], zernike_mode.size))
        for i in range(0, zernike_mode.size):
            z = self._zg.getZernike(zernike_mode[i])
            aa = np.ma.masked_array(z.data, mask=img.mask)
            mat.T[i] = aa.compressed()

        self._mat = mat
        inv = np.linalg.pinv(mat)
        a = np.dot(inv, img.compressed())
        return a, mat

    def _zernikeSurface(self, surface_zernike_coeff_array, ima_mask,
                       mat, index=None):
        '''
        Parameters
        ----------
            surface_zernike_coeff_array: numpy array
                                        vector containing the amplitudes
                                        of the Zernike modes
            ima_mask: numpy array
                    the mask area in which we want to rebuilding
                    the surfaces
            mat: numpy array
                interaction matrix for the masked area
            index: int
                    vector containing the index number of interaction matrix
                    that we want to use
        Returns
        -------
            surf: numpy masked array
                reconstructed surface
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
