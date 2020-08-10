'''
@author: cselmi
'''

import os
import numpy as np
from astropy.io import fits as pyfits
from m4.utils.interface_4D import comm4d
from m4.configuration.config import fold_name
from m4.ground import tracking_number_folder
from m4.utils.zernike_on_m_4 import ZernikeOnM4
from m4.utils.parabola_identification import ParabolIdent

class RotOptAlign():

    def __init__(self, ott):
        """The constructor """
        self._c4d = comm4d()
        self._ott = ott
        self._zOnM4 = ZernikeOnM4()
        self._parab = ParabolIdent()

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return fold_name.ROT_OPT_ALIGN_ROOT_FOLDER


    def acquire_image(self):
        save = tracking_number_folder.TtFolder(RotOptAlign._storageFolder())
        dove, tt = save._createFolderToStoreMeasurements()
        self._ott.angle(0)
        rot_angle = 20
        number_of_image = np.int(360/rot_angle)
        angle_list = []
        self._cube = None

        for k in range(number_of_image):
            star_angle = self._ott.angle()
            self._ott.angle(star_angle + rot_angle)
            angle_list.append(star_angle + rot_angle)
            masked_ima = self._c4d.acq4d(self._ott, 1, show=1)
            name = 'Frame_%04d.fits' %k
            self._saveInterfData(dove, name, masked_ima)
            if self._cube is None:
                self._cube = masked_ima
            else:
                self._cube = np.ma.dstack((self._cube, masked_ima))
        self._saveCube(dove)
        self._saveAngles(angle_list, tt)
        return tt

    def analyser(self, tt):
        cube = self._readCube(tt)
        tip, tilt = self._tipTiltCalculator(cube)
        centro, axs, raggio = self._parab._fitEllipse(tip, tilt)

    def _saveInterfData(self, dove, file_name, image):
        fits_file_name = os.path.join(dove, file_name)
        pyfits.writeto(fits_file_name, image.data)
        pyfits.append(fits_file_name, image.mask.astype(int))

    def _saveAngles(self, angle_list, tt):
        fits_file_name = os.path.join(RotOptAlign._storageFolder(), tt, 'angle_list.txt')
        file = open(fits_file_name, 'w+')
        file.write('Angles used \n')
        for angle in angle_list:
            file.write('%s \n' %angle)
        file.close()

    def _saveCube(self, dove):
        fits_file_name = os.path.join(dove, 'Cube.fits')
        pyfits.writeto(fits_file_name, self._cube.data)
        pyfits.append(fits_file_name, self._cube.mask.astype(int))

    def _readCube(self, tt):
        file_name = os.path.join(RotOptAlign._storageFolder(), tt, 'Cube.fits')
        hduList = pyfits.open(file_name)
        self._cube = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
        return self._cube

    def _tipTiltCalculator(self, cube):
        coef_tilt_list = []
        coef_tip_list = []
        for i in range(cube.shape[2]):
            image = cube[:,:,i]
            coef, mat = self._zOnM4.zernikeFit(image,
                                               np.array([2, 3]))
            coef_tip_list.append(coef[0])
            coef_tilt_list.append(coef[1])
        tip = np.array(coef_tip_list)
        tilt = np.array(coef_tilt_list)
        return tip, tilt

    def _rotationMatrix(self, theta):
        rot_mat = np.zeros((2,2))
        rot_mat[0,0] = np.cos(theta)
        rot_mat[0,1] = - np.sin(theta)
        rot_mat[1,0] = np.sin(theta)
        rot_mat[1,1] = np.cos(theta)
        return rot_mat
        