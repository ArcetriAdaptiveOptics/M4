'''
@author: cs
'''

import numpy as np
from m4.utils.optical_alignment import opt_alignment
from m4.utils.optical_calibration import opt_calibration
from m4.ground import object_from_fits_file_name as obj
from m4.utils.roi import ROI
from m4.utils.zernike_on_m_4 import ZernikeOnM4

class Alignment():
    """
    Class to be used for alignment of the optical tower
    and the deformable mirror

    HOW TO USE IT:
    from m4.alignment import Alignment
    a = Alignment()
    #for the optical tower
    tt = a.ott_calibration(commandAmpVector, nPushPull, maskIndex)
    cmd = a.ott_alignement(tt)
    #for deformable mirror
    tt, zCoefComa, comaSurface = a.m4_calibration(...)
    cmd = a.m4_alignement(zCoefComa)
    """

    def __init__(self):
        """The constructor """
        self._cal = opt_calibration()
        self._roi = ROI()
        self._zOnM4 = ZernikeOnM4()


    def ott_calibration(self, command_amp_vector, n_push_pull, mask_index):
        '''
            Calibration of the optical tower
            args:
                command_amp_vector = vector containing the movement values
                                    of the 5 degrees of freedom
                n_push_pull = number of push pull for each degree of freedom
                mask_index = int (3 for the RM mask)

            returns:
                    tt = tracking number of measurement
        '''
        self._moveRM()
        self._tt = self._cal.measureCalibrationMatrix(0, command_amp_vector,
                                                      n_push_pull)
        int_mat, rec = self._cal.analyzerCalibrationMeasurement(self._tt,
                                                                mask_index)
        return self._tt

    def ott_alignement(self, tt=None):
        """
        Args:
            tt = tracking number of measurement of which you want to use the
                interaction matrix and reconstructor
                None for the last measurement made
        Returns:
                cmd = vector of command to apply to PAR+RM dof
        """
        if tt is None:
            al = opt_alignment(self._tt)
        else:
            al = opt_alignment(tt)
        cmd = al.opt_align()
        self._applyCmd()
        return cmd


    def m4_calibration(self, commandAmpVector_ForParRmAlignement,
                       nPushPull_ForParRmAlignement,
                       maskIndex_ForParRmAlignement,
                       commandAmpVector_ForM4Calibration,
                       nPushPull_ForM4Calibration, maskIndex_ForM4Alignement):
        """
        Calibration of the deformable mirror

        Args:
            commandAmpVector_ForParRmAlignement = amplitude to be applied to par+rm
            nPushPull_ForParRmAlignemen = number of push pull for par+rm dof
            maskIndex_ForParRmAlignement = number of mask index to use
            commandAmpVector_ForM4Calibration = amplitude to be applied to m4
            nPushPull_ForM4Calibration = number of push pull for m4 dof
            maskIndex_ForM4Alignement = number of mask index to use

        Returns:
            zernike_coef_coma = zernike coefficient value for coma calculated
                                by the function _measureComaOnSegmentMask
            coma_surface = reconstructed surface
        """
        self._moveSegmentView()
        tt = self.ott_calibration(commandAmpVector_ForParRmAlignement,
                                  nPushPull_ForParRmAlignement,
                                  maskIndex_ForParRmAlignement)
        cmd = self.ott_alignement(tt)
        self._moveRM()
        zernike_coef_coma, coma_surface = self._measureComaOnSegmentMask()
        self._tt = self._cal.measureCalibrationMatrix(3,
                                                      commandAmpVector_ForM4Calibration,
                                                      nPushPull_ForM4Calibration)
        intMat, rec = self._cal.analyzerCalibrationMeasurement(self._tt,
                                                               maskIndex_ForM4Alignement)
        return self._tt, zernike_coef_coma, coma_surface

    def m4_alignement(self, zernike_coef_coma, tt=None):
        """
        Args:
            zernike_coef_coma = zernike coefficient value for coma
            tt = tracking number of measurement of which you want to use the
                interaction matrix and reconstructor
                None for the last measurement made
        Returns:
                cmd = vector of command to apply to M4 dof
        """
        if tt is None:
            al = opt_alignment(self._tt)
        else:
            al = opt_alignment(tt)
        cmd = al.opt_align(zernike_coef_coma)
        return cmd

    def _moveRM(self):
        pass

    def _moveSegmentView(self):
        pass

    def _applyCmd(self):
        pass

    def _measureComaOnSegmentMask(self):
        ima = obj.readImageFromFitsFileName('Allineamento/20191001_081344/img.fits')
        roi = self._roi.roiGenerator(ima)
        segment_ima = np.ma.masked_array(ima.data, mask=roi[11])

        coef, mat = self._zOnM4.zernikeFit(segment_ima, np.arange(2, 11))
        coma = coef[5]
        coma_surface = self._zOnM4.zernikeSurface(np.array([coma]),
                                                  segment_ima.mask, mat,
                                                  np.array([5]))
        return coma, coma_surface
