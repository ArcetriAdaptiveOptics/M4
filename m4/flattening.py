'''
@author: cs
'''

import os
import logging
import numpy as np
import pyfits
from m4.ground.configuration import Configuration
from m4.ground.tracking_number_folder import TtFolder
from m4.utils.img_redux import TipTiltDetrend


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
        self._command = None
        self._flatteningWf = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                            "Flattening")

    def readVMatrix(self):
        """ Function that returns V matrix (892, 811) for the segment
        """
        root = Configuration.V_MATRIX_FOR_SEGMENT_ROOT_811
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
        self._logger.info('Calculation of the flat command')
        self._an.setDetectorMask(wf.mask | self._an.getMasterMask())
        rec = self._an.getReconstructor()
        wf_masked = np.ma.masked_array(wf.data,
                                       mask=np.ma.mask_or(wf.mask,
                                                          self._an.getMasterMask()))
        amp = -np.dot(rec, wf_masked.compressed())
        v_matrix_cut = self.readVMatrix()
        self._command = np.dot(v_matrix_cut, amp)
        return self._command

    def syntheticWfCreator(self, wf_mask, command):
        """ Returns the synthetic wavefront using the input command
        and the same masked wavefront used to determine the command itself.
        """
        v_matrix_cut = self.readVMatrix()
        v_pinv = np.linalg.pinv(v_matrix_cut)
        amp = np.dot(v_pinv, command)
        sintetic_wf = np.dot(self._an.getInteractionMatrix(), amp)

        mm = np.ma.mask_or(wf_mask, self._an.getMasterMask())
        final_wf_data = np.zeros((Configuration.DIAMETER_IN_PIXEL_FOR_SEGMENT_IMAGES,
                                  Configuration.DIAMETER_IN_PIXEL_FOR_SEGMENT_IMAGES))
        final_wf_data[np.where(mm == False)] = sintetic_wf
        final_wf = np.ma.masked_array(final_wf_data, mask=mm)
        return final_wf

    def flattening(self, offset):
        """ Function for flat command application
        """
        self._logger.info('Application of the flat command')
        #misuro la posizione dello specchio (pos)
        pos = np.zeros(7)

        cmd = pos + self._command
        #la applico e misuro il nuovo wf
        pass



    def save(self):
        store_in_folder = Flattenig._storageFolder()
        save = TtFolder(store_in_folder)
        dove, tt = save._createFolderToStoreMeasurements()
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
