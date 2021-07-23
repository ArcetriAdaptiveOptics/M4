'''
@author: cs
'''

import os
import logging
import numpy as np
from astropy.io import fits as pyfits
from m4.configuration.ott_parameters import OttParameters
from m4.utils.roi import ROI
#from m4.configuration.create_ott import DMirror
from m4.ground.timestamp import Timestamp
from m4.configuration.config import fold_name
from m4.ground import tracking_number_folder


class Flattenig():
    """ Class dealing with the determination and application of the wave front
    flattening command.
    """

    def __init__(self, analyzerIFFunctions):
        """The constructor:
            analyzerIFFunctions = analyzer object to use
        """
        self._logger = logging.getLogger('FLATTENING:')
        self._an = analyzerIFFunctions
        self._who = self._an._who
        self._roi = ROI()
        self._mirror = DMirror()
        self._command = None
        self._flatteningWf = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return fold_name.FLAT_ROOT_FOLD

    def readVMatrix(self):
        """ Function that returns V matrix (892, 811) for the segment
        """
        root = OttParameters.V_MATRIX_FOR_SEGMENT_ROOT_811
        #root = Configuration.V_MATRIX_FOR_SEGMENT_ROOT_892
        hduList = pyfits.open(root)
        v_matrix = hduList[0].data
        v_matrix_cut = v_matrix[:, 0:811]
        return v_matrix_cut

#comando che permette di ottenere la misura del wf dall'interferometro (wf)
#ampr = np.random.randn(25)
#wf = np.dot(self._an._cube, ampr)
    def flatCommand(self, wf):
        """ Returns the command to give to the actuators to level
        the wf considered
        """
        tt = Timestamp.now()
        fitsFileName = os.path.join(Flattenig._storageFolder(), tt)

        self._logger.info('Calculation of the flat command')
        circular_mask = self._roi._circularMaskForSegmentCreator()
        mask_no_edge_actuators = np.ma.mask_or(wf.mask, circular_mask)

        normal_mask = np.ma.mask_or(wf.mask, self._an.getMasterMask())
        super_mask = np.ma.mask_or(mask_no_edge_actuators, self._an.getMasterMask())
        wf_masked = np.ma.masked_array(wf.data, mask=super_mask)
        pyfits.writeto(os.path.join(fitsFileName, 'imgstart.fits'), wf_masked.data)
        pyfits.append(os.path.join(fitsFileName, 'imgstart.fits'), wf_masked.mask.astype(int))

        self._an.setDetectorMask(wf_masked.mask | self._an.getMasterMask())
        rec = self._an.getReconstructor()

        amp = -np.dot(rec, wf_masked.compressed())
        v_matrix_cut = self.readVMatrix()
        self._command = np.dot(v_matrix_cut, amp)
        pyfits.writeto(os.path.join(fitsFileName, 'flatDeltaCommand.fits'), self._command)
        return self._command

    def syntheticWfCreator(self, wf_mask, command):
        """ Returns the synthetic wavefront using the input command
        and the same masked wavefront used to determine the command itself.
        """
        v_matrix_cut = self.readVMatrix()
        v_pinv = np.linalg.pinv(v_matrix_cut)
        amp = np.dot(v_pinv, command)
        sintetic_wf = np.dot(self._an.getInteractionMatrix(), amp)

        circular_mask = self._roi._circularMaskForSegmentCreator()
        mask_no_edge_actuators = np.ma.mask_or(wf_mask, circular_mask)
        normal_mm = np.ma.mask_or(wf_mask, self._an.getMasterMask())
        super_mm = np.ma.mask_or(mask_no_edge_actuators, self._an.getMasterMask())
        final_wf_data = np.zeros((OttParameters.DIAMETER_IN_PIXEL_FOR_SEGMENT_IMAGES,
                                  OttParameters.DIAMETER_IN_PIXEL_FOR_SEGMENT_IMAGES))
        final_wf_data[np.where(super_mm == False)] = sintetic_wf
        final_wf = np.ma.masked_array(final_wf_data, mask=super_mm)
        return final_wf

    def flattening(self, command, seg):
        """ Function for flat command application
        """
        self._logger.info('Application of the flat command')
        #dentro delta command misura già la posizione dello specchio (pos)
        delta_command = self._mirror.mirror_command(command, seg)
        #la applico e misuro il nuovo wf
        return delta_command



    def save(self):
        """ Function to save the info file
        """
        store_in_folder = Flattenig._storageFolder()
        dove, tt = tracking_number_folder.createFolderToStoreMeasurements(store_in_folder)
        fits_file_name = os.path.join(dove, 'info.fits')
        header = pyfits.Header()
        header['WHO'] = self._who
        pyfits.writeto(fits_file_name, self._command, header)
        pyfits.append(fits_file_name, self._flatteningWf.data, header)
        pyfits.append(fits_file_name, self._flatteningWf.mask.astype(int),
                      header)
        return tt


    @staticmethod
    def load(tracking_number):
        """ Function for reload
        """
        theObject = Flattenig(tracking_number)
        store_in_folder = Flattenig._storageFolder()
        folder = os.path.join(store_in_folder, tracking_number)
        info_fits_file_name = os.path.join(folder, 'info.fits')
        header = pyfits.getheader(info_fits_file_name)
        hduList = pyfits.open(info_fits_file_name)
        theObject._who = header['WHO']
        theObject._command = hduList[0].data
        theObject._flattening = np.ma.masked_array(hduList[1].data,
                                                   hduList[2].data.astype(bool))
        return theObject




#### Usato solo qui quindi temporaneamente posizionato qui ###
from m4.configuration.config import path_name
from m4.configuration import ott_parameters

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

