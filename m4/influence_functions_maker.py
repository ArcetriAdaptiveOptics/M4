'''
@author: cs
'''

import os
import copy
import logging
import h5py
import numpy as np
import pyfits
from m4.ground import tracking_number_folder
from m4.ground.configuration import Configuration
from m4.type.commandHistory import CmdHistory
from m4.ground import object_from_fits_file_name


class IFFunctionsMaker(object):

    def __init__(self, device):
        self._device = device
        self._who = self._device._who
        self._nActs = self._device.nActs()
        self._logger= logging.getLogger('IFF_MAKER:')

    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                            "IFFunctions")


    def acq_IFFunctions(self, modes_vector_fits_file_name, n_push_pull,
                        amplitude_fits_file_name, cmd_matrix_fits_file_name,
                        shuffle=None):

        amplitude, modes_vector, cmd_matrix = \
        object_from_fits_file_name.readObjectFitsFileName(amplitude_fits_file_name,
                                                          modes_vector_fits_file_name,
                                                          cmd_matrix_fits_file_name)
        self._nPushPull = n_push_pull

        self._modesVectorTag = modes_vector_fits_file_name
        self._amplitudeTag = amplitude_fits_file_name
        self._cmdMatrixTag = cmd_matrix_fits_file_name
        '''
         arg:
             modesVectorTag= fits file name (modes.fits)
             n_push_pull= numero di push pull consecutivi sull'attuatore
                         (int)
             amplitudeTag= fits file name (amp.fits)
             cmdMatrixTag= fits file name (modalBase.fits)
             shuffle= se non indicato viene creata la matrice della storia
                     dei comandi ordinata

        return:
                tracking number delle misure effettuate
        '''
        store_in_folder = self._storageFolder()
        indexing_input = copy.copy(modes_vector)
        save = tracking_number_folder.TtFolder(store_in_folder)
        dove, tt = save._createFolderToStoreMeasurements()

        diagonal_mat = self._diagonalControll(cmd_matrix)
        if np.count_nonzero(diagonal_mat -
                            np.diag(np.diagonal(diagonal_mat))) == 0:
            print('Measure of zonal IF')
            self._logger.info("Measurement of zonal influence functions for %s. Location: %s",
                       self._who, tt)

        else:
            print('Measure of global IF')
            self._logger.info("Measurement of modal influence functions for %s. Location: %s",
                       self._who, tt)


        cmdH = CmdHistory(self._device)

        if shuffle is None:
            command_history_matrix_to_apply, self._tt_cmdH = \
                    cmdH.tidyCommandHistoryMaker(modes_vector,
                                                 amplitude,
                                                 cmd_matrix,
                                                 n_push_pull)
        else:
            command_history_matrix_to_apply, self._tt_cmdH = \
                    cmdH.shuffleCommandHistoryMaker(modes_vector,
                                                    amplitude,
                                                    cmd_matrix,
                                                    n_push_pull)
        self._indexingList = cmdH.getIndexingList()
        self._saveInfoAsFits(dove)


        self._applyToDM()

        return tt


    def _diagonalControll(self, matrix):
        v = np.zeros((self._nActs, 1))
        reps = matrix.shape[0] - matrix.shape[1]
        vects = np.tile(v, reps)
        new_matrix = np.hstack((matrix, vects))

        return new_matrix


    def _applyToDM(self):
        pass

    def _testIFFunctions_createCube25fromFileFitsMeasure(self):
        fold = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/OIM_25modes.fits'
        hduList = pyfits.open(fold)
        cube_50 = hduList[0].data

        imaList = []
        maskList = []
        for i in range(cube_50.shape[0]):
            if i%2 == 0:
                imaList.append(cube_50[i])
            else:
                maskList.append(cube_50[i])

        cube_25 = None
        zipped = zip(imaList, maskList)
        for ima, mask in zipped:
            immagine = np.ma.masked_array(ima, mask=mask)
            if cube_25 is None:
                cube_25 = immagine
            else:
                cube_25 = np.ma.dstack((cube_25, immagine))
        return cube_25


    def _saveInfoAsFits(self, folder):
        fits_file_name = os.path.join(folder, 'info.fits')
        header = pyfits.Header()
        header['NPUSHPUL'] = self._nPushPull
        header['WHO'] = self._who
        header['TT_CMDH'] = self._tt_cmdH
        header['MODEVECT'] = self._modesVectorTag
        header['CMDMAT'] = self._cmdMatrixTag
        header['AMP'] = self._amplitudeTag
        pyfits.writeto(fits_file_name, self._indexingList, header)

    def _saveInfoAsH5(self, folder):
        fits_file_name = os.path.join(folder, 'info.h5')
        hf = h5py.File(fits_file_name, 'w')
        hf.create_dataset('dataset_1', data=self._indexingList)
        hf.attrs['MODEVECT'] = self._modesVectorTag
        hf.attrs['CMDMAT'] = self._cmdMatrixTag
        hf.attrs['AMP'] = self._amplitudeTag
        hf.attrs['NPUSHPUL'] = self._nPushPull
        hf.attrs['WHO'] = self._who
        hf.attrs['TT_CMDH'] = self._tt_cmdH
        hf.close()


    @staticmethod
    def loadInfoFromFits(folder):
        additional_info_fits_file_name = os.path.join(folder, 'info.fits')
        header = pyfits.getheader(additional_info_fits_file_name)
        hduList = pyfits.open(additional_info_fits_file_name)
        acts_vector_tag = header['MODEVECT']
        cmd_matrix_tag = header['CMDMAT']
        cmd_ampl_tag = header['AMP']
        indexingList = hduList[0].data
        who = header['WHO']
        tt_cmdH = header['TT_CMDH']
        try:
            n_push_pull = header['NPUSHPUL']
        except KeyError:
            n_push_pull = 1

        cmd_ampl, acts_vector, cmd_matrix = \
        object_from_fits_file_name.readObjectFitsFileName(cmd_ampl_tag,
                                                          acts_vector_tag,
                                                          cmd_matrix_tag)
        return (who, tt_cmdH, acts_vector, cmd_matrix, cmd_ampl,
                n_push_pull, indexingList)

    @staticmethod
    def loadInfoFromH5(tt):
        store_in_folder = IFFunctionsMaker._storageFolder()
        folder = os.path.join(store_in_folder, tt)
        file_name = os.path.join(folder, 'info.h5')
        hf = h5py.File(file_name, 'r')
        hf.keys()
        data_1 = hf.get('dataset_1')
        acts_vector_tag = hf.attrs['MODEVECT']
        cmd_matrix_tag = hf.attrs['CMDMAT']
        cmd_ampl_tag = hf.attrs['AMP']
        indexingList = np.array(data_1)
        n_push_pull = hf.attrs['NPUSHPUL']
        who = hf.attrs['WHO']
        tt_cmdH = hf.attrs['TT_CMDH']

        cmd_ampl, acts_vector, cmd_matrix = \
        object_from_fits_file_name.readObjectFitsFileName(cmd_ampl_tag,
                                                          acts_vector_tag,
                                                          cmd_matrix_tag)
        return (who, tt_cmdH, acts_vector, cmd_matrix, cmd_ampl,
                n_push_pull, indexingList)
