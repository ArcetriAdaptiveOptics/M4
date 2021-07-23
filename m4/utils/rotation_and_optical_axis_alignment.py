'''
Authors
  - C. Selmi:  written in September 2020
'''

import os
import time
import numpy as np
import logging
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.utils import tip_tilt_interf_fit
from m4.configuration.config import fold_name
from m4.ground import tracking_number_folder
from m4.ground import zernike
from m4.utils.parabola_identification import ParabolIdent
from m4.noise_functions import Noise
from m4.configuration.ott_parameters import OpcUaParameters

class RotOptAlign():
    """
    Class for alignment data acquisition and analysis::

        fom m4.configuration import start
        ott = start.create_ott()
        from m4.utils.rotation_and_optical_axis_alignment import RotOptAlign
        ro = RotOptAlign(ott)
        tt = ro.image_acquisition(start_point, end_point, n_points)
        centro, axs, raggio = ro.data_analyzer(tt)
    """
    def __init__(self, ott, interf):
        """The constructor """
        self._logger = logging.getLogger('ROTOPTALIGN')
        self._c4d = interf
        self._ott = ott
        self._parab = ParabolIdent()
        self._n = Noise()

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return fold_name.ROT_OPT_ALIGN_ROOT_FOLDER

    def _checkAngle(self, angle):
        if angle <= OpcUaParameters.min_angle or angle >= OpcUaParameters.max_angle:
            raise OSError(' The required angle is incorrect: %d' % angle)
        else:
            pass

    def image_acquisition(self, start_point, end_point, n_points):
        """
        Parameters
        ----------
                start_point: int
                            value of start angle
                end_point: int
                            value of end angle
                n_points:int
                        number of images desired

        Returns
        -------
                tt: string
                    tracking number of measurements
        """
        self._logger.info('Images acquisition')
        dove, tt = tracking_number_folder.createFolderToStoreMeasurements(RotOptAlign._storageFolder())

        self._checkAngle(start_point)
        self._checkAngle(end_point)

        total_angle = np.abs(end_point - start_point)/1.
        self._ott.angle(start_point)
        rot_angle = total_angle/n_points
        number_of_image = np.int(total_angle/rot_angle)
        angle_list = []
        self._cube = None

        start_angle = self._ott.angle()
        if end_point < start_point:
            direction = -1
        elif end_point > start_point:
            direction = 1

        for k in range(number_of_image+1):
            #start_angle = self.ott.angle()
            self._ott.angle(start_angle + k*rot_angle*direction)
            angle_list.append(start_angle + k*rot_angle*direction)
            time.sleep(5)
            masked_ima = self._c4d.acq4d(1, self._ott)
            name = 'Frame_%04d.fits' %k
            self._c4d.save_phasemap(dove, name, masked_ima)
            if self._cube is None:
                self._cube = masked_ima
            else:
                self._cube = np.ma.dstack((self._cube, masked_ima))
        self._saveCube(dove)
        self._saveAngles(angle_list, tt)
        return tt

    def data_analyzer(self, tt):
        """
        Parameters
        ----------
                tt: string
                    tracking number of measurements

        Returns
        -------
            centro: numpy array
                    coordinates of the center
            axs: numpy array
                major and minor axis coming from the fit of the ellipse
            raggio: int
                    radius of the parabola circumference
        """
        self._logger.info('Images analysis')
        cube = self._readCube(tt)
        tip, tilt = self._tipTiltCalculator(cube)
        self._plot(tip, tilt)
        centro, axs, raggio = self._parab._fitEllipse(tip, tilt)
        return centro, axs, raggio

    def _plot(self, tip, tilt):
        plt.figure(figsize=(7,7))
        plt.plot(tip, tilt, '-o')

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
            coef, mat = zernike.zernikeFit(image,
                                            np.array([2, 3]))
            coef_tip_list.append(-coef[0]/(633e-9) *np.sqrt(2))
            coef_tilt_list.append(coef[1]/(633e-9) *np.sqrt(2))
        tip = np.array(coef_tip_list)
        tilt = np.array(coef_tilt_list)
        #devo diventare angoli
        #mettendo la maschera
        #fare il plot di tip, tilt (plot isometrico:stesse dimensioni x y)
        return tip, tilt

    def singleImage(self, masked_ima=None):
        """
        On single image returns the coefficients of tip and tilt both in Zernike of Noll
        and in units of the interferometer

        Other Parameters
        ---------------
                masked_ima: numpy masked array
                    if it is not passed acquires the interferometer

        Returns
        -------
            tip_tilt: numpy array
                    Zernike Noll coefficients
            cc: numpy array
                    Interferometers values
        """
        if masked_ima is None:
            masked_ima = self._c4d.acq4d(1, self._ott)
        else:
            masked_ima = masked_ima
        coef, mat = zernike.zernikeFit(masked_ima,
                                       np.array([2, 3]))
        tip = -coef[0]/(633e-9) *np.sqrt(2)
        tilt = coef[1]/(633e-9) *np.sqrt(2)

        cc = tip_tilt_interf_fit.fit(masked_ima)
        return np.array([tip, tilt]), cc

    def _rotationMatrix(self, theta):
        rot_mat = np.zeros((2,2))
        rot_mat[0,0] = np.cos(theta)
        rot_mat[0,1] = - np.sin(theta)
        rot_mat[1,0] = np.sin(theta)
        rot_mat[1,1] = np.cos(theta)
        return rot_mat
