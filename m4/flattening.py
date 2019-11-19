'''
@author: cs
'''

import os
import numpy as np
import pyfits
from m4.ground.configuration import Configuration
from m4.ground.tracking_number_folder import TtFolder
from m4.utils.img_redux import TipTiltDetrend


class Flattenig():

    def __init__(self, analyzerIFFunctions):
        self._an = analyzerIFFunctions
        self._who = self._an._who
        self._command = None
        self._flatteningWf = None

    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                            "Flattening")


    def flatCommand(self, wf):
        #comando che permette di ottenere la misura del wf dall'interferometro (wf)
        #ampr = np.random.randn(25)
        #wf = np.dot(self._an._cube, ampr)

        self._an.setDetectorMask(wf.mask | self._an.getMasterMask())
        rec = self._an.getReconstructor()
        wf_masked = np.ma.masked_array(wf.data,
                                       mask=np.ma.mask_or(wf.mask,
                                                          self._an.getMasterMask()))
        amp = -np.dot(rec, wf_masked.compressed())
        return amp


    def flattening(self, offset):
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
        theObject._who= header['WHO']
        theObject._command= hduList[0].data
        theObject._flattening= np.ma.masked_array(hduList[1].data,
                                                  hduList[2].data.astype(bool))
        return theObject
