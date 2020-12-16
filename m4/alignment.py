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
from m4.utils.zernike_on_m_4 import ZernikeOnM4
from m4.utils.interface_4D import comm4d

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
        self._zOnM4 = ZernikeOnM4()
        self._ott = ott
        self._c4d = comm4d()


    def ott_calibration(self, command_amp_vector, n_push_pull, mask_index):
        '''Calibration of the optical tower

        Parameters
        ----------
                command_amp_vector: numpy array
                                  vector containing the movement values
                                  of the 5 degrees of freedom
                n_push_pull: int
                            number of push pull for each degree of freedom
                mask_index: int
                            3 for the RM mask (non ruotate 2)

        Returns
        -------
                tt: string
                    tracking number of measurements made
        '''
        self._tt = self._cal.measureCalibrationMatrix(self._ott, 0, command_amp_vector,
                                                      n_push_pull)
        int_mat, rec = self._cal.analyzerCalibrationMeasurement(self._tt,
                                                                mask_index)
        return self._tt

    def ott_alignment(self, n_images, move, intMatModesVector=None, tt=None):
        """
        Parameters
        ----------
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
        par_cmd, rm_cmd, dove = al.opt_align(self._ott, n_images, intMatModesVector)
        if move == 1:
            self._write_par(par_cmd)
            self._write_rm(rm_cmd)
            image = self._c4d.acq4d(self._ott, n_images)
            name = 'FinalImage.fits'
            self._c4d.save_phasemap(dove, name, image)
        return par_cmd, rm_cmd


    def m4_calibration(self, commandAmpVector_ForM4Calibration,
                       nPushPull_ForM4Calibration, maskIndex_ForM4Alignement):
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
        zernike_coef_coma, coma_surface = self._measureComaOnSegmentMask()
        print(zernike_coef_coma)
        self._tt = self._cal.measureCalibrationMatrix(self._ott, 3,
                                                      commandAmpVector_ForM4Calibration,
                                                      nPushPull_ForM4Calibration)
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

    def _moveRM(self, rslide):
        self._logger.debug('rslide = %f', rslide)
        self._ott.rslide(rslide)

    def _moveSegmentView(self, slide, angle):
        self._logger.debug('slide = %f angle = %f', slide, angle)
        self._ott.slide(slide)
        self._ott.angle(angle)

    def _write_par(self, par_cmd):
        pos_par = self._ott.parab()
        self._ott.parab(pos_par + par_cmd)

    def _write_rm(self, rm_cmd):
        pos_rm = self._ott.refflat()
        self._ott.refflat(pos_rm + rm_cmd)
        #p,m = self._c4d.acq4d(self._ott, 1, show=1)

    def _write_m4(self, m4_cmd):
        pos_m4 = self._ott.m4()
        self._ott.m4(pos_m4 + m4_cmd)
        #p,m = self._c4d.acq4d(self._ott, 1, show=1)

    def _measureComaOnSegmentMask(self):
        #ima = obj.readImageFromFitsFileName('Allineamento/20191001_081344/img.fits')
        ima = self._c4d.acq4d(self._ott, 1, show=1)
        roi = self._roi.roiGenerator(ima)
        segment_ima = np.ma.masked_array(ima.data, mask=roi[5])

        coef, mat = self._zOnM4.zernikeFit(segment_ima, np.arange(2, 11))
        coma = coef[5]
        coma_surface = self._zOnM4.zernikeSurface(np.array([coma]),
                                                  segment_ima.mask, mat,
                                                  np.array([5]))
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

# fare un fits e via
    def _readZcoef(self, tt):
        dove = os.path.join(self._cal._storageFolder(), tt)
        fits_file_name = os.path.join(dove, 'z_coma.fits')
        hduList = pyfits.open(fits_file_name)
        header = pyfits.getheader(fits_file_name)
        z_coma = header['COMA']
        return z_coma
