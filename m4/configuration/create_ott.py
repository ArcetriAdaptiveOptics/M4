'''
Authors
  - C. Selmi: written in 2020
'''

import logging
import os
import numpy as np
from m4.configuration import config as conf, ott_parameters
from m4.configuration.config import path_name
from m4.configuration.ott_parameters import OttParameters, OpcUaParameters
from m4.utils.roi import ROI
from m4.ground.opc_ua_controller import OpcUaController


class OTT():

    def __init__(self, parabola_slider, reference_mirror_slider, angle_rotator):
        """The constructor """
        self._logger = logging.getLogger('OTT:')
        self._r = ROI()
        self._opcUa = OpcUaController()
        self._parabola_slider = parabola_slider
        self._reference_mirror_slider = reference_mirror_slider
        self._angle_rotator = angle_rotator

        self.par_start_position = np.zeros(6)
        self.m4_start_position = np.zeros(6)
        self.refflat_start_position = np.zeros(6)


# Elements position
    def slide(self, par_trans=None):
        ''' Function to set the parabola translation (range: -0.9 m +0.9 m)

        Other Parameters
        ----------
        par_trans: int, optional [mm]
                If par_trans is not set it's equal to zero

        Returns
        -------
            par_trans: int [mm]
                    parabola translation
        '''
        if par_trans is None:
            return self._parabola_slider.getPosition()
        else:
            return self._parabola_slider.setPosition(par_trans)

    def rslide(self, ref_flat=None):
        '''  Function to set the reference flat mirror (range: -0.05 m to 0.4 m)

        Other Parameters
        ----------
        ref_flat: int, optional [mm]
                If ref_flat is not set it's equal to zero

        Returns
        -------
        ref_flat: int [mm]
                reference flat mirror position
        '''
        if ref_flat is None:
            return self._reference_mirror_slider.getPosition()
        else:
            return self._reference_mirror_slider.setPosition(ref_flat)

    def angle(self, rot_ring_angle=None):
        ''' Function to set the rotating ring angle (range: 0 to 360)

        Other Parameters
        ----------
            rot_ring_angle: int, optional [deg]
                If rot_ring_angle is not set it's equal to zero

        Returns
        -------
            rot_ring_angle: int [deg]
                            rot_ring_angle
        '''
        if rot_ring_angle is None:
            return self._angle_rotator.getPosition()
        else:
            return self._angle_rotator.setPosition(rot_ring_angle)


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
        # if type(start_position) is np.ndarray:
        # if start_position.size == 6:
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
                    # print(start_position[j])
                self._opcUa.move_object(OpcUaParameters.PAR_KIN)
                self._opcUa.wait_for_stop(OpcUaParameters.PAR_KIN)
                self.par_start_position = self._readParPosition()
        # else:
            # raise OSError('Incorrect length of the vector')
        # else:
            # raise OSError('Data is not a numpy array')
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
        # if type(start_position) is np.ndarray:
        # if start_position.size == 6:
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
                    # print(start_position[j])
                self._opcUa.move_object(OpcUaParameters.RM_KIN)
                self._opcUa.wait_for_stop(OpcUaParameters.RM_KIN)
                self.refflat_start_position = self._readRMPosition()
        # else:
            # raise OSError('Incorrect length of the vector')
        # else:
            # raise OSError('Data is not a numpy array')
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
        # if type(start_position) is np.ndarray:
        # if start_position.size == 6:
        if conf.simulated == 1:
            if start_position is None:
                self.m4_start_position = self.m4_start_position
            else:
                self.m4_start_position = start_position
            self._logger.debug(self.m4_start_position)
        else:
            print('Sw to be developed')
            self.m4_start_position = self.m4_start_position
        # else:
            # raise OSError('Incorrect length of the vector')
        # else:
            # raise OSError('Data is not a numpy array')
        return self.m4_start_position

    def temperature(self):
        if conf.simulated == 1:
            temp_vector = np.zeros(OpcUaParameters.num_PT_sensor)
        else:
            temp_vector = self._opcUa.get_temperature_vector()
        return temp_vector



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
            if command_input.shape[0] == OttParameters.N_ACTS_TOT:
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
        # forza = self._mirror._ff * delta_command
        return delta_command

    def _measurePosition(self):
        # dall'opc ua va letta la posizione degli attuatori
        pos = np.zeros(OttParameters.N_ACTS_TOT) + 7
        return pos
