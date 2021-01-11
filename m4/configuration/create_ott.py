'''
Authors
  - C. Selmi: written in 2020
'''

import logging
import os
import numpy as np
from astropy.io import fits as pyfits
from m4.configuration import config as conf, ott_parameters
from m4.configuration.config import path_name
from m4.configuration.ott_parameters import OttParameters, Interferometer, OpcUaParameters
from m4.ground import read_data
from m4.utils.roi import ROI
from m4.ground import zernike, geo
from m4.ground.opc_ua_controller import OpcUaController

class OTT():

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('OTT:')
        self._opcUa = OpcUaController()
        self._r = ROI()
#        self._zg = ZernikeGenerator(2*OttParameters.parab_radius*OttParameters.pscale)
        self._slide = 0
        self._rslide = 0
        self._angle = 0
        self.par_start_position = np.zeros(6)
        self.m4_start_position = np.zeros(6)
        self.refflat_start_position = np.zeros(6)
        self.m4offset = 0.
        self.offset = 0.
        self.smap = np.zeros((Interferometer.N_PIXEL[0], Interferometer.N_PIXEL[1]))
        self.rmap = np.zeros(((2*OttParameters.rflat_radius*OttParameters.pscale).astype(int),
                              (2*OttParameters.rflat_radius*OttParameters.pscale).astype(int)))
        self.m4pupil = read_data.readFits_data(os.path.join(conf.path_name.MIRROR_FOLDER,
                                                        conf.mirror_conf,
                                                        'm4_mech_pupil-bin2.fits'))
        self.m4ima = self.m4pupil * 0.
        self.mask = read_data.readFits_data(os.path.join(conf.path_name.MIRROR_FOLDER,
                                                     conf.mirror_conf,
                                                     'ott_mask.fits'))
        self.parmask = np.ma.make_mask(read_data.readFits_data(
                                        os.path.join(conf.path_name.OPTICAL_FOLDER,
                                        conf.optical_conf, 'ottmask.fits')))
#         self.segmask1 = np.ma.make_mask(obj.readFits_object(
#                                         os.path.join(conf.path_name.MIRROR_FOLDER,
#                                         conf.mirror_conf, 'py-sect4-mask.fits')))
#         self.ifmat = obj.readFits_object(os.path.join(conf.path_name.MIRROR_FOLDER,
#                                                     conf.mirror_conf,
#                                                    'if_sect4_rot-bin2.fits'))
#         self.vmat = obj.readFits_object(os.path.join(conf.path_name.MIRROR_FOLDER,
#                                                    conf.mirror_conf, 'ff_v_matrix.fits'))
        self.zmat = read_data.readFits_data(os.path.join(conf.path_name.OPTICAL_FOLDER,
                                                     conf.optical_conf,
                                                     'Zmat.fits'))

# Elements position
    def slide(self, par_trans=None):
        ''' Function to set the parabola translation (range: -0.9 m +0.9 m)

        Other Parameters
        ----------
        par_trans: int, optional
                If par_trans is not set it's equal to zero

        Returns
        -------
            par_trans: int
                    parabola translation
        '''
        self._logger.debug('About SLIDE')
        if conf.simulated == 1:
            if par_trans is None:
                self._slide = self._slide
            else:
                self._slide = par_trans
            self._logger.debug('Position = %f', self._slide)
        else:
            if par_trans is None:
                self._slide = self._opcUa.get_position(OpcUaParameters.ST)
            else:
                self._checkSlide(par_trans)
                self._slide = self._opcUa.set_target_position(OpcUaParameters.ST, par_trans)
                self._opcUa.move_object(OpcUaParameters.ST)
                self._opcUa.wait_for_stop(OpcUaParameters.ST)
                self._slide = self._opcUa.get_position(OpcUaParameters.ST)
        return self._slide

    def _checkSlide(self, slide):
        if slide <= OpcUaParameters.min_slide or slide >= OpcUaParameters.max_slide:
            raise OSError(' The required parabola position is incorrect: %d' % slide)
        else:
            pass

    def rslide(self, ref_flat=None):
        '''  Function to set the reference flat mirror (range: -0.05 m to 0.4 m)

        Other Parameters
        ----------
        ref_flat: int, optional
                If ref_flat is not set it's equal to zero

        Returns
        -------
        ref_flat: int
                reference flat mirror position
        '''
        self._logger.debug('About RSLIDE')
        if conf.simulated == 1:
            if ref_flat is None:
                self._rslide = self._rslide
            else:
                self._rslide = ref_flat
            self._logger.debug('Position = %f', self._rslide)
        else:
            if ref_flat is None:
                self._rslide = self._opcUa.get_position(OpcUaParameters.CAR)
            else:
                self._checkRslide(ref_flat)
                self._rslide = self._opcUa.set_target_position(OpcUaParameters.CAR, ref_flat)
                self._opcUa.move_object(OpcUaParameters.CAR)
                self._opcUa.wait_for_stop(OpcUaParameters.CAR)
                self._rslide = self._opcUa.get_position(OpcUaParameters.CAR)
        return self._rslide

    def _checkRslide(self, r_slide):
        if r_slide <= OpcUaParameters.min_r_slide or r_slide >= OpcUaParameters.max_r_slide:
            raise OSError(' The required reference flat position is incorrect: %d' % r_slide)
        else:
            pass

    def angle(self, rot_ring_angle=None):
        ''' Function to set the rotating ring angle (range: 0 to 360)

        Other Parameters
        ----------
            rot_ring_angle: int, optional
                If rot_ring_angle is not set it's equal to zero

        Returns
        -------
            rot_ring_angle: int
                            rot_ring_angle
        '''
        self._logger.debug('About ANGLE')
        if conf.simulated == 1:
            if rot_ring_angle is None:
                self._angle = self._angle
            else:
                self._angle = rot_ring_angle
            self._logger.debug('Position = %f', self._angle)
        else:
            if rot_ring_angle is None:
                self._angle = self._opcUa.get_position(OpcUaParameters.RA)
            else:
                self._checkAngle(rot_ring_angle)
                self._opcUa.set_target_position(OpcUaParameters.RA, rot_ring_angle)
                self._opcUa.move_object(OpcUaParameters.RA)
                self._opcUa.wait_for_stop(OpcUaParameters.RA)
                self._angle = self._opcUa.get_position(OpcUaParameters.RA)
        return self._angle

    def _checkAngle(self, angle):
        if angle <= OpcUaParameters.min_angle or angle >= OpcUaParameters.max_angle:
            raise OSError(' The required angle is incorrect: %d' % angle)
        else:
            pass

# Elements alignment
    def parab(self, start_position=None):
        '''Function to set the start position of the parable

        Other Parameters
        ----------
        start_position: numpy array, optional
                        vector of six position,
                        If start_position is not set it's equal to zero vector

        Returns
        -------
            start_position: numpy array
                        start position of the parable

        '''
        self._logger.debug('About PARAB')
        #if type(start_position) is np.ndarray:
        #if start_position.size == 6:
        if conf.simulated == 1:
            if start_position is None:
                self.par_start_position = self.par_start_position
            else:
                self.par_start_position = start_position
            self._logger.debug(self.par_start_position)
        else:
            if start_position is None:
                self.par_start_position = self._readParPosition()
            else:
                n_opc = np.array([OpcUaParameters.PAR_PISTON,
                                  OpcUaParameters.PAR_TIP,
                                  OpcUaParameters.PAR_TILT])
                for i in range(OttParameters.PARABOLA_DOF.size):
                    j = OttParameters.PARABOLA_DOF[i]
                    self._opcUa.set_target_position(n_opc[i], start_position[j])
                    #print(start_position[j])
                self._opcUa.move_object(OpcUaParameters.PAR_KIN)
                self._opcUa.wait_for_stop(OpcUaParameters.PAR_KIN)
                self.par_start_position = self._readParPosition()
        #else:
            #raise OSError('Incorrect length of the vector')
        #else:
            #raise OSError('Data is not a numpy array')
        return self.par_start_position

    def _readParPosition(self):
        piston = self._opcUa.get_position(OpcUaParameters.PAR_PISTON)
        tip = self._opcUa.get_position(OpcUaParameters.PAR_TIP)
        tilt = self._opcUa.get_position(OpcUaParameters.PAR_TILT)
        return np.array([0, 0, piston, tip, tilt, 0])

    def refflat(self, start_position=None):
        '''Function to set the start position of the reference flat

        Other Parameters
        ----------
        start_position: numpy array, optional
                        vector of six position,
                        If start_position is not set it's equal to zero vector

        Returns
        -------
            start_position: numpy array
                        start position of the reference flat
        '''
        self._logger.debug('About REFFLAT')
        #if type(start_position) is np.ndarray:
        #if start_position.size == 6:
        if conf.simulated == 1:
            if start_position is None:
                self.refflat_start_position = self.refflat_start_position
            else:
                self.refflat_start_position = start_position
            self._logger.debug(self.refflat_start_position)
        else:
            if start_position is None:
                self.refflat_start_position = self._readRMPosition()
            else:
                n_opc = np.array([OpcUaParameters.RM_PISTON,
                                  OpcUaParameters.RM_TIP,
                                  OpcUaParameters.RM_TILT])
                for i in range(OttParameters.RM_DOF_PISTON.size):
                    j = OttParameters.RM_DOF_PISTON[i]
                    self._opcUa.set_target_position(n_opc[i], start_position[j])
                    #print(start_position[j])
                self._opcUa.move_object(OpcUaParameters.RM_KIN)
                self._opcUa.wait_for_stop(OpcUaParameters.RM_KIN)
                self.refflat_start_position = self._readRMPosition()
        #else:
            #raise OSError('Incorrect length of the vector')
        #else:
            #raise OSError('Data is not a numpy array')
        return self.refflat_start_position

    def _readRMPosition(self):
        piston = self._opcUa.get_position(OpcUaParameters.RM_PISTON)
        tip = self._opcUa.get_position(OpcUaParameters.RM_TIP)
        tilt = self._opcUa.get_position(OpcUaParameters.RM_TILT)
        return np.array([0, 0, piston, tip, tilt, 0])

    def m4(self, start_position=None):
        '''Function to set the start position of the deformable mirror

        Other Parameters
        ----------
        start_position: numpy array, optional
                        vector of six position,
                        If start_position is not set it's equal to zero vector

        Returns
        -------
            start_position: numpy array
                        start position of the deformable mirror
        '''
        self._logger.debug('About M4')
        #if type(start_position) is np.ndarray:
        #if start_position.size == 6:
        if conf.simulated == 1:
            if start_position is None:
                self.m4_start_position = self.m4_start_position
            else:
                self.m4_start_position = start_position
            self._logger.debug(self.m4_start_position)
        else:
            print('Sw to be developed')
            self.m4_start_position = self.m4_start_position
        #else:
            #raise OSError('Incorrect length of the vector')
        #else:
            #raise OSError('Data is not a numpy array')
        return self.m4_start_position

### Sensitivity matrices
    def _readMatFromTxt(self, file_name):
        ''' Function to read matrix of 11 Zernike x 6 displacements,
        m RMS, per 1 m displacement - or 1 radiant rotation

        Parameters
        ----------
        file_name: string
                    matrix file path

        Returns
        -------
                mat: numpy array [11,6]
                    matrix from txt file
        '''
        file = open(file_name, 'r')
        triplets = file.read().split()
        x = np.array(triplets)
        mat = x.reshape(11, 6)
        return mat.astype(float)

    def zmx_parpos2z(self):
        '''
        Returns
        -------
                mat: numpy array [11,6]
                    matrix parable positions to zernike
        '''
        file_name = os.path.join(conf.path_name.OPTICAL_FOLDER,
                                 conf.optical_conf, 'PAR_pos2z.txt')
        #file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ZST_PAR_pos2z.txt'
        mat = self._readMatFromTxt(file_name)
        return mat

    def zmx_refflatpos2z(self):
        '''
        Returns
        -------
                mat: numpy array [11,6]
                    matrix reference flat positions to zernike
        '''
        file_name = os.path.join(conf.path_name.OPTICAL_FOLDER,
                                 conf.optical_conf, 'M4_pos2z.txt')
        #file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ZST_FM_pos2z.txt'
        mat = self._readMatFromTxt(file_name)
        return mat

    def zmx_m4pos2z(self):
        '''
        Returns
        -------
                mat: numpy array [11,6]
                    matrix deformable mirror positions to zernike
        '''
        file_name = os.path.join(conf.path_name.OPTICAL_FOLDER,
                                 conf.optical_conf, 'M4_pos2z.txt')
        #file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ZST_M4_pos2z.txt'
        mat = self._readMatFromTxt(file_name)
        return mat
# Zmat
    def create_zmat(self, file_name):
        '''
        Returns
        -------
            zmat: numpy array

        '''
        #file_name = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/OTT/ott_mask.fits'
        hduList = pyfits.open(file_name)
        final_mask = np.invert(hduList[0].data.astype(bool))

        prova = np.ma.masked_array(np.ones(((2*OttParameters.parab_radius*OttParameters.pscale).astype(int),
                                            (2*OttParameters.parab_radius*OttParameters.pscale).astype(int))),
                                   mask=final_mask)
        zernike_mode = np.arange(10)+1
        mm = np.where(final_mask == False)
        x, y, r, xx, yy = geo.qpupil(final_mask)
        zmat = zernike.getZernike(xx[mm], yy[mm], zernike_mode)
        return zmat

class DMirror():
    def __init__(self):
        """The constructor """
        curr_conffolder = os.path.join(path_name.CONFIGURATION_ROOT_FOLDER,
                                       ott_parameters.tnconf_mirror)
#         self.vmat = read_data.readFits_data(os.path.join(curr_conffolder, 'vmat.fits'))
#         self.ff = read_data.readFits_data(os.path.join(curr_conffolder, 'ff_matrix.fits'))

        self.m4od = OttParameters.m4od
        self.m4optod = OttParameters.m4optod
        self.m4id = OttParameters.m4id
        self._activeSegment = 0

    def mirror_command(self, command, seg=None):
        command_input = np.copy(command)
        pos = self._measurePosition()
        if seg is None:
            if command_input.shape[0]==OttParameters.N_ACTS_TOT:
                command = command_input
            elif command_input.shape[0] < OttParameters.N_ACTS_TOT:
                cmd = np.zeros(OttParameters.N_ACTS_TOT)
                for j in range(command_input.shape[0]):
                    act = j + (OttParameters.N_ACT_SEG * self._activeSegment)
                    cmd[act] = command_input[j]
                command = cmd
            delta_command = pos + command
        else:
            command_list = []
            for i in range(seg.shape[0]):
                cmd = np.zeros(OttParameters.N_ACTS_TOT)
                for j in range(OttParameters.N_ACT_SEG):
                    act = j + (OttParameters.N_ACT_SEG * seg[i])
                    k = i * OttParameters.N_ACT_SEG
                    cmd[act] = command_input[k]
                    command = cmd
                command_list.append(command)
            command = np.zeros(OttParameters.N_ACTS_TOT)
            for cmd in command_list:
                command = command + cmd
            delta_command = command
        #forza = self._mirror._ff * delta_command
        return delta_command

    def _measurePosition(self):
        #dall'opc ua va letta la posizione degli attuatori
        pos = np.zeros(OttParameters.N_ACTS_TOT)+7
        return pos
