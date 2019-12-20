'''
@author: cselmi
'''

import os
import logging
import numpy as np
from m4.ground import tracking_number_folder
from m4.ground.configuration import Configuration

class SPL():

    def __init__(self, filter_obj, camera_obj):
        """The constructor """
        self._logger = logging.getLogger('SPL_CONTROLLER:')
        self._filter = filter_obj
        self._camera = camera_obj

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                            "SPL")

    def measurement(self):
        ''' Deve contenere tutta la procedura di acquisizione ed analisi dei dati
        '''
        self._exptime, self._acq4d, self._an = self._setParameters(0.7, 1, 1)


    def _setParameters(self, exptime, acq4d, an):
        self._exptime = exptime
        self._camera.ExposureTimeAbs = 1.5 * exptime * 1e6
        self._camera.Timeout = 30
        self._acq4d = acq4d
        self._an = an
        return self._exptime, self._acq4d, self._an

    def acquire(self):
        ''' Acquisizione dati 
        '''
        store_in_folder = self._storageFolder()
        save = tracking_number_folder.TtFolder(store_in_folder)
        dove, tt = save._createFolderToStoreMeasurements()



    def _preProcessing(self, image):
        ''' Prepara i dati acquisiti per essere analizzati
        '''
        counts, bin_edges = np.histogram(image)
        thr = 5 * bin_edges(counts == max(counts))
        idx = np.where(image < thr)
        image(idx) = 0


    def analyzer(self):
        ''' Analizza i dati e li confronta con quelli sintetici
        '''

