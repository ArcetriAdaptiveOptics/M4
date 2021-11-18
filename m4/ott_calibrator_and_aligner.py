'''
Authors
  - C. Selmi: written in 2019
              modified in 2021
'''

import os
import logging
from m4.utils.optical_alignment import OpticalAlignment
from m4.utils.optical_calibration import OpticalCalibration
from m4.utils.roi import ROI


class OttCalibAndAlign():
    """
    Class to be used for alignment of the optical tower
    and the deformable mirror

    HOW TO USE IT::

        from m4.alignment import OttCalibAndAlign
        from m4.configuration import start
        ott, interf = start.create_ott(conf='.../youConf.yaml')
        ac = OttCalibAndAlign(ott, interf)
        #for PAR+RM
        tt_calib = a.par_and_rm_calibrator(commandAmpVector, nPushPull, maskIndex)
        par_cmd, rm_cmd = a.par_and_rm_aligner(tt_calib)
    """

    def __init__(self, ott, interf):
        """The constructor """
        self._logger = logging.getLogger('ALIGNMENT:')
        self._ott = ott
        self._interf = interf
        self._cal = OpticalCalibration(ott, interf)
        self._tt = None
        self._roi = ROI()

    def par_and_rm_calibrator(self, n_frames, command_amp_vector, n_push_pull):
        '''Calibration of the optical tower

        Parameters
        ----------
        command_amp_vector: numpy array
                          vector containing the movement values
                          of the 5 degrees of freedom
        n_push_pull: int
                    number of push pull for each degree of freedom
        n_frames: int
                number of frame for 4D measurement

        Returns
        -------
        tt: string
            tracking number of measurements made
        '''
        self._tt = self._cal.measureAndAnalysisCalibrationMatrix(0, command_amp_vector,
                                                                 n_push_pull, n_frames)
        return self._tt

    def par_and_rm_aligner(self, move, tt_cal, n_images,
                      zernike_to_be_corrected=None, dof_command_id=None):
        """
        Parameters
        ----------
            n_images: int
                number of interferometers frames
            move: boolean
                True to move the tower
                other to show commands
            tt: string, None
                tracking number of measurement of which you want to use the
                interaction matrix and reconstructor

        Other Parameters
        ----------
        zernike_to_be_corrected: numpy array
                        None is equal to np.array([0,1,2,3,4,5])
                        for tip, tilt, fuoco, coma, coma
        dof_command_id: numpy array
                array containing the number of degrees of freedom to be commanded

        Returns
        -------
                par_cmd: numpy array
                    vector of command to apply to PAR dof
                rm_cmd: numpy array
                    vector of command to apply to RM dof
        """
        aliner = OpticalAlignment(tt_cal, self._ott, self._interf)
        par_cmd, rm_cmd, dove = aliner.opt_aligner(n_images,
                                                   zernike_to_be_corrected,
                                                   dof_command_id)
        if move is True:
            pos_par = self._ott.parabola.getPosition()
            self._ott.parabola.setPosition(pos_par + par_cmd)
            pos_rm = self._ott.referenceMirror.getPosition()
            self._ott.referenceMirror.setPosition(pos_rm + rm_cmd)
        image = self._interf.acquire_phasemap(n_images)
        name = 'FinalImage.fits'
        all_final_coef, final_coef_selected = aliner.getZernikeWhitAlignerObjectOptions(image)
        self._alignmentLog(aliner, all_final_coef, dof_command_id, move)
        self._interf.save_phasemap(dove, name, image)
        return par_cmd, rm_cmd, dove

    def _alignmentLog(self, aligner, total_coef, dof_command_id, move):
        fits_file_name = os.path.join(aligner._storageFolder(), 'AlignmentLog.txt')
        file = open(fits_file_name, 'a+')
        for i in range(total_coef.size):
            file.write('%9.3e ' % total_coef[i])
        file.write('\n')
        if move == 0:
            dof_command_id = -1
        file.write('%s \n ************\n' % dof_command_id)
        file.close()


### M4 calibrator and aligner in cartellaBella.m4.toImplement.ott_calibrator_and_aligner ###
    def m4_calibrator(self):
        pass
    
    def m4_aligner(self):
        pass
