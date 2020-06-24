'''
Autors
  - C. Selmi:  written in 2019
'''

import numpy as np
from m4.utils.optical_alignment import opt_alignment
from m4.utils.optical_calibration import opt_calibration
from m4.ground import object_from_fits_file_name as obj
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
        tt, zCoefComa, comaSurface = a.m4_calibration(...)
        cmd = a.m4_alignement(zCoefComa)
    """

    def __init__(self, ott):
        """The constructor """
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
        self._moveSegmentView(0.75, 90.)
        self._moveRM(0.6)
        self._tt = self._cal.measureCalibrationMatrix(self._ott, 0, command_amp_vector,
                                                      n_push_pull)
        int_mat, rec = self._cal.analyzerCalibrationMeasurement(self._tt,
                                                                mask_index)
        return self._tt

    def ott_alignement(self, tt=None):
        """
        Parameters
        ----------
            tt: string, None
                tracking number of measurement of which you want to use the
                interaction matrix and reconstructor
                None for the last measurement made
        Returns
        -------
                cmd: numpy array
                    vector of command to apply to PAR+RM dof
        """
        if tt is None:
            al = opt_alignment(self._tt)
        else:
            al = opt_alignment(tt)
        par_cmd, rm_cmd = al.opt_align(self._ott)
        self._applyCmd(par_cmd, rm_cmd)
        return par_cmd, rm_cmd


    def m4_calibration(self, commandAmpVector_ForM4Calibration,
                       nPushPull_ForM4Calibration, maskIndex_ForM4Alignement):
        """ Calibration of the deformable mirror

        Parameters
        ----------
            commandAmpVector_ForParRmAlignement: numpy array
                                            amplitude to be applied to par+rm
            nPushPull_ForParRmAlignemen: int
                                        number of push pull for par+rm dof
            maskIndex_ForParRmAlignement: int
                                        number of mask index to use
                                        (reference mirror mask)
            commandAmpVector_ForM4Calibration: numpy array
                                            amplitude to be applied to m4
            nPushPull_ForM4Calibration: int
                                        number of push pull for m4 dof
            maskIndex_ForM4Alignement: int
                                        number of mask index to use
                                        (segment mask)

        Returns
        -------
            zernike_coef_coma: int
                                zernike coefficient value for coma calculated
                                by the function _measureComaOnSegmentMask
            coma_surface: numpy array
                            reconstructed surface
        """
#         self._moveSegmentView(0.75, 90.)
#         self._moveRM(0.6)
#         tt = self.ott_calibration(commandAmpVector_ForParRmAlignement,
#                                   nPushPull_ForParRmAlignement,
#                                   maskIndex_ForParRmAlignement)
#         par_cmd, rm_cmd = self.ott_alignement(tt)
        zernike_coef_coma, coma_surface = self._measureComaOnSegmentMask()
        self._tt = self._cal.measureCalibrationMatrix(self._ott, 3,
                                                      commandAmpVector_ForM4Calibration,
                                                      nPushPull_ForM4Calibration)
        intMat, rec = self._cal.analyzerCalibrationMeasurement(self._tt,
                                                               maskIndex_ForM4Alignement)
        return self._tt, zernike_coef_coma, coma_surface

    def m4_alignement(self, zernike_coef_coma, tt=None):
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
                cmd: numpy array
                    vector of command to apply to M4 dof
        """
        self._moveRM(0.)
        if tt is None:
            al = opt_alignment(self._tt)
        else:
            al = opt_alignment(tt)
        cmd = al.opt_align(self._ott, zernike_coef_coma)
        return cmd

    def _moveRM(self, rslide):
        self._ott.rslide(rslide)

    def _moveSegmentView(self, slide, angle):
        self._ott.slide(slide)
        self._ott.angle(angle)

    def _applyCmd(self, par_cmd, rm_cmd):
        pos_par = self._ott.parab()
        self._ott.parab(pos_par + par_cmd)
        pos_rm = self._ott.refflat()
        self._ott.refflat(pos_rm + rm_cmd)

    def _measureComaOnSegmentMask(self):
        #ima = obj.readImageFromFitsFileName('Allineamento/20191001_081344/img.fits')
        p, m = self._c4d.acq4d(self._ott, 1, show=1)
        ima = np.ma.masked_array(p.T, mask=np.invert(m.astype(bool)).T)
        roi = self._roi.roiGenerator(ima)
        segment_ima = np.ma.masked_array(ima.data, mask=roi[5])

        coef, mat = self._zOnM4.zernikeFit(segment_ima, np.arange(2, 11))
        coma = coef[5]
        coma_surface = self._zOnM4.zernikeSurface(np.array([coma]),
                                                  segment_ima.mask, mat,
                                                  np.array([5]))
        return coma, coma_surface
