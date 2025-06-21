"""
Authors
  - C. Selmi:  written in September 2020
"""

import os
import time
import logging
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.configuration import folders as fold_name
from m4.ground import tracking_number_folder
from m4.ground import zernike
from m4.utils.parabola_identification import ParabolIdent


class RotOptAlign:
    """
    Class for alignment data acquisition and analysis::

        from m4.configuration import start
        conf = '..../MyConfiguration.yaml'
        ott, interf = start.create_ott(conf)
        from m4.utils.rotation_and_optical_axis_alignment import RotOptAlign
        ro = RotOptAlign(ott, interf)
        tt = ro.image_acquisition(start_point, end_point, n_points)
        centro, axs, raggio = ro.data_analyzer()
    """

    def __init__(self, ott, interf):
        """The constructor"""
        self._logger = logging.getLogger("ROTOPTALIGN")
        self._interf = interf
        self._ott = ott
        self._parab = ParabolIdent()
        self._cube = None
        self._angles = None
        self.tt = None

    @staticmethod
    def _storageFolder():
        """Creates the path where to save measurement data"""
        return fold_name.ROT_OPT_ALIGN_ROOT_FOLDER

    def image_acquisition(self, start_point, end_point, n_points):
        """
        Parameters
        ----------
        start_point: int [deg]
                    absolute position in deg for the value of start angle
        end_point: int [deg]
                    absolute position in deg for the value of end angle
        n_points:int
                number of images desired

        Returns
        -------
        tt: string
            tracking number of measurements
        """
        self._logger.info("Images acquisition")
        dove, self.tt = tracking_number_folder.createFolderToStoreMeasurements(
            RotOptAlign._storageFolder()
        )

        total_angle = np.abs(end_point - start_point) / 1.0
        self._ott.angleRotator.setPosition(start_point)
        rot_angle = total_angle / n_points
        number_of_image = int(total_angle / rot_angle)
        angle_list = []
        self._cube = None

        start_angle = self._ott.angleRotator.getPosition()
        if end_point < start_point:
            direction = -1
        elif end_point > start_point:
            direction = 1

        for k in range(number_of_image + 1):
            # start_angle = self.ott.angle()
            self._ott.angleRotator.setPosition(start_angle + k * rot_angle * direction)
            angle_list.append(start_angle + k * rot_angle * direction)
            time.sleep(5)
            masked_ima = self._interf.acquire_phasemap(1)
            name = "Frame_%04d.fits" % k
            self._interf.save_phasemap(dove, name, masked_ima)
            if self._cube is None:
                self._cube = masked_ima
            else:
                self._cube = np.ma.dstack((self._cube, masked_ima))
        self._saveCube(dove)
        self._angles = angle_list
        self._saveAngles(self._angles, self.tt)
        return self.tt

    def data_analyzer(self):
        """
        Returns
        -------
            centro: numpy array
                    coordinates of the center
            axs: numpy array
                major and minor axis coming from the fit of the ellipse
            raggio: int
                    radius of the parabola circumference
        """
        self._logger.info("Images analysis")
        cube = self.getCube()
        tip, tilt = self._tipTiltCalculator(cube)
        self._plot(tip, tilt)
        centro, axs, raggio = self._parab._fitEllipse(tip, tilt)
        return centro, axs, raggio

    def _plot(self, tip, tilt):
        plt.figure(figsize=(7, 7))
        plt.plot(tip, tilt, "-o")

    def _saveAngles(self, angle_list, tt):
        fits_file_name = os.path.join(
            RotOptAlign._storageFolder(), tt, "angle_list.txt"
        )
        file = open(fits_file_name, "w+")
        file.write("Angles used \n")
        for angle in angle_list:
            file.write("%s \n" % angle)
        file.close()
        fits_file_name = os.path.join(RotOptAlign._storageFolder(), tt, "theta.fits")
        pyfits.writeto(fits_file_name, np.array(angle_list))

    def _saveCube(self, dove):
        fits_file_name = os.path.join(dove, "Cube.fits")
        pyfits.writeto(fits_file_name, self._cube.data)
        pyfits.append(fits_file_name, self._cube.mask.astype(int))

    def _readTheta(self, tt):
        file_name = os.path.join(RotOptAlign._storageFolder(), tt, "theta.fits")
        hduList = pyfits.open(file_name)
        self._theta = hduList[0].data
        return self._theta

    def _readCube(self, tt):
        file_name = os.path.join(RotOptAlign._storageFolder(), tt, "Cube.fits")
        hduList = pyfits.open(file_name)
        self._cube = np.ma.masked_array(
            hduList[0].data, mask=hduList[1].data.astype(bool)
        )
        return self._cube

    def getCube(self):
        """
        Returns
        -------
        cube: numpy masked array
            image acquired during the measurements
        """
        if self._cube is None:
            self._cube = self._readCube(self.tt)
        return self._cube

    def getAngles(self):
        """
        Returns
        -------
        angles: numpy array
            values of angles [deg] used during the measuremnts
        """
        if self._angles is None:
            self._angles = self._readTheta(self.tt)
        return self._angles

    def _tipTiltCalculator(self, cube):
        coef_tilt_list = []
        coef_tip_list = []
        for i in range(cube.shape[2]):
            image = cube[:, :, i]
            coef, mat = zernike.zernikeFit(image, np.array([2, 3]))
            coef_tip_list.append(-coef[0] / (633e-9) * np.sqrt(2))
            coef_tilt_list.append(coef[1] / (633e-9) * np.sqrt(2))
        tip = np.array(coef_tip_list)
        tilt = np.array(coef_tilt_list)
        # devo diventare angoli
        # mettendo la maschera
        # fare il plot di tip, tilt (plot isometrico:stesse dimensioni x y)
        return tip, tilt

    def _rotationMatrix(self, theta):
        rot_mat = np.zeros((2, 2))
        rot_mat[0, 0] = np.cos(theta)
        rot_mat[0, 1] = -np.sin(theta)
        rot_mat[1, 0] = np.sin(theta)
        rot_mat[1, 1] = np.cos(theta)
        return rot_mat

    @staticmethod
    def reloadROObject(tt):
        """Function used to create the rotation and optical alignment object
        from tracking number

        Parameters
        ----------
        tt: string
            tracking number of measuremnts
        """
        theObject = RotOptAlign("nulla", "niente")
        theObject.tt = tt
        theObject._cube = theObject.getCube()
        theObject._angles = theObject.getAngles()
        return theObject
