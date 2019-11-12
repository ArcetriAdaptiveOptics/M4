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

    def __init__(self):
        self._cal = opt_calibration()
        self._roi = ROI()
        self._zOnM4 = ZernikeOnM4()


    def ott_calibration(self, command_amp_vector, n_push_pull, mask_index):
        '''
            arg:
                command_amp_vector= vettore contenente i valori dei movimenti
                                    dei 5 gradi di libertà
                n_push_pull= numero di push pull per ogni grado di libertà
                mask_index= int (3 per la maschera dell'RM)
        '''
        self._moveRM()
        self._tt = self._cal.measureCalibrationMatrix(0, command_amp_vector,
                                                      n_push_pull)
        int_mat, rec = self._cal.analyzerCalibrationMeasurement(self._tt,
                                                                mask_index)
        return self._tt

    def ott_alignement(self, tt=None):
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
