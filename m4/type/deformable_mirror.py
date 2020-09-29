'''
@author: cs
'''
import numpy as np
import os
from m4.configuration.ott_parameters import *
from m4.ground import object_from_fits_file_name as obj


class m4():
    """ Class for m4 definition"""
    def __init__(self):
        """The constructor """
        self._nActsTot = OttParameters.N_ACTS_TOT
        self._who = 'All segments'

    def nActs(self):
        return self._nActsTot


class segment(m4):
    """ Class for segment definition"""
    def __init__(self, segment_index):
        """The constructor """
        super().__init__()

        if segment_index < 6:
            self._segmentIndex = segment_index
        else:
            raise OSError('Segment number %s doesnt exists' % segment_index)

        self._nActSeg = OttParameters.N_ACT_SEG
        self._nSeg = OttParameters.N_SEG
        self._who = 'Segment number %s' % segment_index
        self._angleInDegrees = OttParameters.REFERENCE_ANGLE_DEGREES
        self._angleInRad = OttParameters.REFERENCE_ANGLE_RAD
        self._distance = OttParameters.SEGMENT_DISTANCE_FROM_CENTRE

    def nActs(self):
        return self._nActSeg


class Mirror():
    def __init__(self):
        """The constructor """
        self._nActsTot = OttParameters.N_ACTS_TOT
        self._nActSeg = OttParameters.N_ACT_SEG
        self._nSeg = OttParameters.N_SEG
        self._activeSegment = 0
        #self._ff = obj.readFits_object(os.path.join(Configuration.curr_conffolder,'ff_matrix.fits'))

    def mirror_command(self, command, seg=None):
        command_input = np.copy(command)
        pos = self._measurePosition()
        if seg is None:
            if command_input.shape[0]==self._nActsTot:
                command = command_input
            elif command_input.shape[0] < self._nActsTot:
                cmd = np.zeros(self._nActsTot)
                for j in range(command_input.shape[0]):
                    act = j + (self._nActSeg * self._activeSegment)
                    cmd[act] = command_input[j]
                command = cmd
            delta_command = pos + command
        else:
            command_list = []
            for i in range(seg.shape[0]):
                cmd = np.zeros(self._nActsTot)
                for j in range(self._nActSeg):
                    act = j + (self._nActSeg * seg[i])
                    k = i * self._nActSeg
                    cmd[act] = command_input[k]
                    command = cmd
                command_list.append(command)
            command = np.zeros(self._nActsTot)
            for cmd in command_list:
                command = command + cmd
            delta_command = command
        #forza = self._mirror._ff * delta_command
        return delta_command

    def _measurePosition(self):
        pos = np.zeros(self._nActsTot)
        return pos