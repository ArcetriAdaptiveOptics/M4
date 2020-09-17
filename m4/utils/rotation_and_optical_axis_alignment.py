'''
@author: cselmi
'''

import os
import time
import numpy as np
import logging
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.utils.interface_4D import comm4d
from m4.configuration.config import fold_name
from m4.ground import tracking_number_folder
from m4.utils.zernike_on_m_4 import ZernikeOnM4
from m4.utils.parabola_identification import ParabolIdent
from m4.noise_functions import Noise

class RotOptAlign():

    def __init__(self, ott):
        """The constructor """
        self._logger = logging.getLogger('ROTOPTALIGN')
        self._c4d = comm4d()
        self._ott = ott
        self._zOnM4 = ZernikeOnM4()
        self._parab = ParabolIdent()
        self._n = Noise()

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return fold_name.ROT_OPT_ALIGN_ROOT_FOLDER


    def acquire_image(self, n_points, start_point, direction):
        self._logger.info('Images acquisition')
        save = tracking_number_folder.TtFolder(RotOptAlign._storageFolder())
        dove, tt = save._createFolderToStoreMeasurements()
        self._ott.angle(start_point)
        #n_points = 24
        rot_angle = 350/n_points
        number_of_image = np.int(350/rot_angle)
        angle_list = []
        self._cube = None

        start_angle = self._ott.angle()
        #direction = -1
        for k in range(number_of_image+1):
            #start_angle = self._ott.angle()
            self._ott.angle(start_angle + k*rot_angle*direction)
            angle_list.append(start_angle + k*rot_angle*direction)
            time.sleep(5)
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

    def analyzer(self, tt):
        self._logger.info('Images analysis')
        cube = self._readCube(tt)
        tip, tilt = self._tipTiltCalculator(cube)
        self._plot(tip, tilt)
        centro, axs, raggio = self._parab._fitEllipse(tip, tilt)
        return centro, axs, raggio

    def _plot(self, tip, tilt):
        plt.figure(figsize=(7,7))
        plt.plot(tip*1e6, tilt*1e6, '-o')

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
        fits_file_name = os.path.join(RotOptAlign._storageFolder(), tt, 'theta.fits')
        pyfits.writeto(fits_file_name, np.array(angle_list))

    def _saveCube(self, dove):
        fits_file_name = os.path.join(dove, 'Cube.fits')
        pyfits.writeto(fits_file_name, self._cube.data)
        pyfits.append(fits_file_name, self._cube.mask.astype(int))

    def _readTheta(self, tt):
        file_name = os.path.join(RotOptAlign._storageFolder(), tt, 'theta.fits')
        hduList = pyfits.open(file_name)
        self._theta = hduList[0].data
        return self._theta

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
            image_ex = self._n._imageExtender(image)
            coef, mat = self._zOnM4.zernikeFit(image_ex,
                                               np.array([2, 3]))
            coef_tip_list.append(coef[0])
            coef_tilt_list.append(coef[1])
        tip = np.array(coef_tip_list)
        tilt = np.array(coef_tilt_list)
        #devo diventare angoli
        #mettendo la maschera
        #fare il plot di tip, tilt (plot isometrico:stesse dimensioni x y)
        return tip, tilt

    def _rotationMatrix(self, theta):
        rot_mat = np.zeros((2,2))
        rot_mat[0,0] = np.cos(theta)
        rot_mat[0,1] = - np.sin(theta)
        rot_mat[1,0] = np.sin(theta)
        rot_mat[1,1] = np.cos(theta)
        return rot_mat

    def tipTiltRotation(self, theta, tip, tilt):
        new_tt_mat = np.zeros((26,2))
        aa = np.dstack((tip,tilt))
        bb = aa[0]
        for i in range(theta.size):
            rot_mat = self._rotationMatrix(theta[i])
            new_tt = np.dot(bb[0], rot_mat)
            new_tt_mat[i,:] = new_tt
        return new_tt_mat
            

    def tipTiltCorrector(self, tip_or_tilt):
        rot_coef = self._tipOrTiltRotation(tip_or_tilt)
        if tip_or_tilt==0:
            posToBeCorrected = 4
        elif tip_or_tilt==1:
            posToBeCorrected = 3
        self._applyCorrection(rot_coef, posToBeCorrected)
        return self._ott.m4()

    def _tipOrTiltRotation(self, tip_or_tilt):
        '''
        Parameters
        ----------
        tip_or_tilt: int
                    0 for tip
                    1 for tilt
        Returns
        -------
        rot_coef: int
        '''
        theta = self._ott.angle()
        masked_ima = self._c4d.acq4d(self._ott, 1, show=1)
        coef, mat = self._zOnM4.zernikeFit(masked_ima,
                                            np.array([2, 3]))
        if tip_or_tilt==0:
            rot_mat = self._rotationMatrix(theta)
        elif tip_or_tilt==1:
            rot_mat = self._rotationMatrix(-theta)
        rot_coef = np.dot(coef, rot_mat)
        return rot_coef[tip_or_tilt]

    def _applyCorrection(self, rot_coef, posToBeCorrected):
        start_position = self._ott.m4()
        new_position = start_position
        new_position[posToBeCorrected] = start_position[posToBeCorrected] + rot_coef
        self._ott.m4(new_position)
