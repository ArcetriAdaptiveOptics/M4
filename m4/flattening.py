'''
@author: cs
'''

import os
import logging
import numpy as np
from astropy.io import fits as pyfits
from m4.configuration.ott_parameters import OttParameters, M4Parameters
from m4.utils import roi
#from m4.configuration.create_ott import DMirror
from m4.ground.timestamp import Timestamp
from m4.configuration import config_folder_names as fold_name
from m4.ground import tracking_number_folder


class Flattenig():
    """ Class dealing with the determination and application of the wave front
    flattening command.
    """

    def __init__(self, analyzerIFFunctions, deformableMirror):
        """The constructor:
            analyzerIFFunctions = analyzer object to use
        """
        self._logger = logging.getLogger('FLATTENING:')
        self._an = analyzerIFFunctions
        self._who = self._an._who
        self._roi = roi
        self._mirror = deformableMirror
        self._command = None
        self._flatteningWf = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return fold_name.FLAT_ROOT_FOLD

    def readVMatrix(self):
        """ Function that returns V matrix (892, 811) for the segment
        """
        root = M4Parameters.V_MATRIX_FOR_SEGMENT_ROOT_811
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
