'''
Authors
  - C. Selmi:  written in 2019
'''

import os
import logging
import numpy as np
from astropy.io import fits as pyfits
from m4.utils.optical_alignment import opt_alignment
from m4.utils.optical_calibration import opt_calibration
from m4.utils.roi import ROI
from m4.ground import zernike
from m4.ground.interface_4D import comm4d

class Alignment():
    """
    Class to be used for alignment of the optical tower
    and the deformable mirror

    HOW TO USE IT::

        from m4.alignment import Alignment
        from m4.configuration import start
        ott = start.create_ott()
        a = Alignment(ott)
        #for the optical tower
        tt = a.ott_calibration(commandAmpVector, nPushPull, maskIndex)
        par_cmd, rm_cmd = a.ott_alignement(tt)
        #for deformable mirror
        tt, zCoefComa, comaSurface = a.m4_calibration(commandAmpVector, nPushPull, maskIndex)
        cmd = a.m4_alignement(zCoefComa, tt)
    """

    def __init__(self, ott):
        """The constructor """
        self._logger = logging.getLogger('ALIGNMENT:')
        self._cal = opt_calibration()
        self._roi = ROI()
        self._ott = ott
        self._c4d = comm4d()


    def ott_calibration(self, n_frames, command_amp_vector, n_push_pull, mask_index):
        '''Calibration of the optical tower

        Parameters
        ----------
                command_amp_vector: numpy array
                                  vector containing the movement values
                                  of the 5 degrees of freedom
                n_push_pull: int
                            number of push pull for each degree of freedom
                mask_index: int
                            3 for the simulatore's RM mask (non ruotate 2)
                            0 for standard mask in real ott

        Returns
        -------
                tt: string
                    tracking number of measurements made
        '''
        self._tt = self._cal.measureCalibrationMatrix(self._ott, 0, command_amp_vector,
                                                      n_push_pull, n_frames)
        int_mat, rec = self._cal.analyzerCalibrationMeasurement(self._tt,
                                                                mask_index)
        return self._tt

    def ott_alignment(self, n_images, move, intMatModesVector=None, commandId=None,
                     tt=None):
        """
        Parameters
        ----------
            n_images: int
                number of interferometers frames
            move: int
                1 to move the tower
                other to show commands
        Other Parameters
        ----------
        intMatModesVecor: numpy array
                        None is equal to np.array([0,1,2,3,4,5])
                        for tip, tilt, fuoco, coma, coma
        commandId: numpy array
                array containing the number of degrees of freedom to be commanded
        tt: string, None
                tracking number of measurement of which you want to use the
                interaction matrix and reconstructor
                None for the last measurement made
        Returns
        -------
                par_cmd: numpy array
                    vector of command to apply to PAR dof
                rm_cmd: numpy array
                    vector of command to apply to RM dof
        """
        if tt is None:
            al = opt_alignment(self._tt)
        else:
            al = opt_alignment(tt)
        par_cmd, rm_cmd, dove = al.opt_align(self._ott, n_images, intMatModesVector, commandId)
        if move == 1:
            pos_par = self._ott.parab()
            self._ott.parab(pos_par + par_cmd)
            pos_rm = self._ott.refflat()
            self._ott.refflat(pos_rm + rm_cmd)
            image = self._c4d.acq4d(self._ott, n_images)
            name = 'FinalImage.fits'
            self._c4d.save_phasemap(dove, name, image)
        return par_cmd, rm_cmd


    def m4_calibration(self, n_frames, commandAmpVector_ForM4Calibration,
                       nPushPull_ForM4Calibration, maskIndex_ForM4Alignement,
                       nFrames):
        """ Calibration of the deformable mirror

        Parameters
        ----------
            commandAmpVector_ForM4Calibration: numpy array
                                            amplitude to be applied to m4
            nPushPull_ForM4Calibration: int
                                        number of push pull for m4 dof
            maskIndex_ForM4Alignement: int
                                        number of mask index to use
                                        rm out = 3, rm in = 5
                                        (segment mask)

        Returns
        -------
            zernike_coef_coma: int
                                zernike coefficient value for coma calculated
                                by the function _measureComaOnSegmentMask
            coma_surface: numpy array
                            reconstructed surface
        """
        zernike_coef_coma, coma_surface = self._measureComaOnSegmentMask(nFrames)
        print(zernike_coef_coma)
        self._tt = self._cal.measureCalibrationMatrix(self._ott, 3,
                                                      commandAmpVector_ForM4Calibration,
                                                      nPushPull_ForM4Calibration, n_frames)
        self._saveZcoef(zernike_coef_coma, coma_surface)
        intMat, rec = self._cal.analyzerCalibrationMeasurement(self._tt,
                                                               maskIndex_ForM4Alignement)
        return self._tt, zernike_coef_coma, coma_surface

    def m4_alignment(self, zernike_coef_coma, tt=None):
        """
        Parameters
        ----------
            zernike_coef_coma: int
                                zernike coefficient value for coma
            tt: string, None
                tracking number of measurement of which you want to use the
                interaction matrix and reconstructor
                None for the last measurement made
        Returns
        -------
                m4_cmd: numpy array
                    vector of command to apply to M4 dof
        """
        #self._moveRM(0.)
        if tt is None:
            al = opt_alignment(self._tt)
        else:
            al = opt_alignment(tt)
        m4_cmd = al.opt_align(self._ott, zernike_coef_coma)
        #self._applyM4Command(m4_cmd)
        return m4_cmd


    def _measureComaOnSegmentMask(self, nFrames):
        #ima = obj.readImageFromFitsFileName('Allineamento/20191001_081344/img.fits')
        ima = self._c4d.acq4d(self._ott, nFrames)
        roi = self._roi.roiGenerator(ima)
        segment_ima = np.ma.masked_array(ima.data, mask=roi[5])

        coef, mat = zernike.zernikeFit(segment_ima, np.arange(10)+1)
        coma = coef[6]
        coma_surface = zernike.zernikeSurface(segment_ima, coma, mat)
        return coma, coma_surface

    def _saveZcoef(self, zernike_coef_coma, coma_surface):
        dove = os.path.join(self._cal._storageFolder(), self._tt)
        fits_file_name = os.path.join(dove, 'z_coma.txt')
        file = open(fits_file_name, 'w+')
        file.write('%4e' %zernike_coef_coma)
        file.close()
        fits_file_name = os.path.join(dove, 'z_coma.fits')
        header = pyfits.Header()
        header['COMA'] = zernike_coef_coma
        pyfits.writeto(fits_file_name, coma_surface.data, header)
        pyfits.append(fits_file_name, coma_surface.mask.astype(int), header)

    def _readZcoef(self, tt):
        dove = os.path.join(self._cal._storageFolder(), tt)
        fits_file_name = os.path.join(dove, 'z_coma.fits')
        hduList = pyfits.open(fits_file_name)
        header = pyfits.getheader(fits_file_name)
        z_coma = header['COMA']
        return z_coma
