'''
@author: cs
'''
import os
import copy
import logging
from astropy.io import fits as pyfits
import h5py
import numpy as np
from m4.ground import tracking_number_folder
from m4.ground.configuration import Configuration


class CmdHistory():
    '''
    Class to create and manage the command history matrix

    HOW TO USE IT:
    from m4.type.commandHistory import CmdHistory
    cmdH = CmdHistory()
    '''

    def __init__(self, device):
        """The constructor """
        self._device = device
        self._logger = logging.getLogger('CMD_HISTORY:')
        self._who = self._device._who 
        self._nActs = self._device.nActs()
        self._modeVector = None
        self._nPushPull = None
        self._cmdMatrix = None
        self._cmdHToApply = None
        self._ampVect = None

    def getCommandHistory(self):
        return self._cmdHToApply

    def getIndexingList(self):
        return self._indexingList

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save data"""
        return os.path.join(Configuration.OPD_DATA_FOLDER,
                            "CommandHistory")


    def tidyCommandHistoryMaker(self, mode_vector, amp_vector,
                                cmd_matrix, temp_nPushPull):
        '''
         arg:
            modesVector = Mode or actuator index vector
                            to apply (numpy.array([]))
            ampVector = amplitude mode vector (numpy.array([]))
            cmdMatrix = mode command matrix
                         (nActs x nModes)
                         diagonal matrix in case of zonal commands
            temp_nPushPull = vector of push pull on the actuator
                             (es. np.array([1, -1, 1]))
         return:
             matrixToApply = tidy command history
                             (nAct, nModes x nPushPusll x 2)
             tt = tracking number
        '''
        self._ampVect = amp_vector
        cmd_history = self._tidyCmdHistory(mode_vector, temp_nPushPull, cmd_matrix)
        aa = np.arange(self._cmdHistory.shape[1])
        bb = np.tile(amp_vector, temp_nPushPull.shape[0])
        zipped = zip(aa, bb)
        matrix_to_apply = self._cmdHistoryToApply(zipped)
        self._cmdHToApply = matrix_to_apply
        tt = self.saveAsFits()
        self._logger.info('Creation of the ordered commandHistoryMatrix %s', tt)
        print(tt)

        return matrix_to_apply, tt


    def shuffleCommandHistoryMaker(self, mode_vector, amp_vector,
                                   cmd_matrix, temp_nPushPull):
        '''
         arg:
            modesVector = Mode or actuator index vector
                            to apply (numpy.array([]))
            ampVector = amplitude mode vector (numpy.array([]))
            cmdMatrix = mode command matrix
                         (nActs x nModes)
                         diagonal matrix in case of zonal commands
            nPushPull = number of consecutive push pull on the actuator
                         (int)
         return:
             matrixToApply = shuffle command history
                             (nAct, nModes x nPushPusll x 2)
             tt = tracking number
        '''
        self._ampVect = amp_vector
        cmd_history, indexing_list = self._shuffleCmdHistory(
            mode_vector, temp_nPushPull, cmd_matrix)
        zipped = self._zippedAmplitude(amp_vector)
        matrix_to_apply = self._cmdHistoryToApply(zipped)
        self._cmdHToApply = matrix_to_apply
        tt = self.saveAsFits()
        self._logger.info('Creation of the shuffle commandHistoryMatrix %s', tt)
        print(tt)

        return matrix_to_apply, tt

    def _shuffleCmdHistory(self, mode_vector, temp_nPushPull, cmd_matrix):
        self._modeVector = copy.copy(mode_vector)
        self._nPushPull = temp_nPushPull.shape[0]
        self._template = temp_nPushPull
        self._cmdMatrix = cmd_matrix

        n_frame = mode_vector.size * temp_nPushPull.shape[0]
        matrix_to_apply = np.zeros((self._nActs, n_frame))

        indexingList = []
        for j in range(temp_nPushPull.shape[0]):
            np.random.shuffle(mode_vector)
            indexingList.append(list(mode_vector))

            cmdList = []
            for i in mode_vector:
                cmd = cmd_matrix[:, i]
                cmdList.append(cmd)

            for i in range(len(cmdList)):
                k = len(cmdList)*j + i
                matrix_to_apply.T[k] = cmdList[i]

        self._cmdHistory = matrix_to_apply
        self._indexingList = np.array(indexingList)

        return self._cmdHistory, self._indexingList


    def _tidyCmdHistory(self, mode_vector, temp_nPushPull, cmd_matrix):
        self._modeVector = copy.copy(mode_vector)
        self._nPushPull = temp_nPushPull.shape[0]
        self._template = temp_nPushPull
        self._cmdMatrix = cmd_matrix
        indList = []
        for i in range(temp_nPushPull.shape[0]):
            indList.append(mode_vector)
        self._indexingList = np.array(indList)
        #self._indexingList= np.tile(mode_vector, temp_nPushPull)

        n_frame = mode_vector.size * temp_nPushPull.shape[0]
        matrix_to_apply = np.zeros((self._nActs, n_frame))

        cmdList = []
        for i in mode_vector:
            cmd = cmd_matrix[:, i]
            cmdList.append(cmd)

        for j in range(temp_nPushPull.shape[0]):
            for i in range(len(cmdList)):
                k = len(cmdList)*j + i
                matrix_to_apply.T[k] = cmdList[i]

        self._cmdHistory = matrix_to_apply

        return self._cmdHistory


    def _amplitudeReorganization(self, indexing_input, indexing_list,
                                 amplitude, temp_nPushPull):
        where = []
        for i in indexing_input:
            for j in range(temp_nPushPull.shape[0]):
                a = np.where(indexing_list[j] == i)
                where.append(a)


        where = np.array(where)
        vect = np.zeros(amplitude.shape[0]*temp_nPushPull.shape[0])

        for i in range(amplitude.shape[0]):
            for k in range(temp_nPushPull.shape[0]):
                p = temp_nPushPull.shape[0] * i + k
                indvect = where[p][0][0]+ indexing_input.shape[0] * k
                vect[indvect] = amplitude[i]

        return vect


    def _zippedAmplitude(self, amplitude):
        aa = np.arange(self._cmdHistory.shape[1])
        reorganized_amplitude = \
                self._amplitudeReorganization(self._modeVector,
                                              self._indexingList,
                                              amplitude,
                                              self._template)
        #zipped= np.dstack((aa, reorganized_amplitude))
        zipped = zip(aa, reorganized_amplitude)
        return zipped

    def _cmdHistoryToApply(self, zipped):
        matrix_with_amp = self._cmdHistory
        for i, amp in zipped:
            matrix_with_amp.T[i] = matrix_with_amp.T[i] * amp

        vec_push_pull = np.array((1, -1))
        matrix_to_apply = np.zeros((self._nActs,
                                    self._cmdHistory.shape[1] *
                                    vec_push_pull.shape[0]))

        for i in range(self._cmdHistory.shape[1]):
            j = 2*i
            matrix_to_apply.T[j] = matrix_with_amp.T[i]* vec_push_pull[0]
            matrix_to_apply.T[j+1] = matrix_with_amp.T[i]* vec_push_pull[1]

        return matrix_to_apply


    def saveAsFits(self):
        """ Save the data in fits format """
        store_in_folder = CmdHistory._storageFolder()
        save = tracking_number_folder.TtFolder(store_in_folder)
        dove, tt = save._createFolderToStoreMeasurements()

        fits_file_name = os.path.join(dove, 'info.fits')
        header = pyfits.Header()
        header['NPUSHPUL'] = self._nPushPull
        header['WHO'] = self._who
        pyfits.writeto(fits_file_name, self._modeVector, header)
        pyfits.append(fits_file_name, self._indexingList, header)
        pyfits.append(fits_file_name, self._cmdMatrix, header)
        pyfits.append(fits_file_name, self._cmdHToApply, header)
        pyfits.append(fits_file_name, self._ampVect, header)
        return tt

    def saveAsH5(self):
        """ Save the data in h5 format """
        store_in_folder = CmdHistory._storageFolder()
        save = tracking_number_folder.TtFolder(store_in_folder)
        dove, tt = save._createFolderToStoreMeasurements()

        fits_file_name = os.path.join(dove, 'info.h5')
        hf = h5py.File(fits_file_name, 'w')
        hf.create_dataset('dataset_1', data=self._modeVector)
        hf.create_dataset('dataset_2', data=self._indexingList)
        hf.create_dataset('dataset_3', data=self._cmdHToApply)
        hf.create_dataset('dataset_4', data=self._ampVect)
        hf.attrs['NPUSHPUL'] = self._nPushPull
        hf.attrs['WHO'] = self._who
        hf.close()
        return tt

    @staticmethod
    def loadFromFits(device, tt):
        """ Creates the object from the info.fits file located in tt

            Args:
                tt = tracking number

            Returns:
                    theObject = CmdHistory class object
        """
        theObject = CmdHistory(device)
        theObject._tt = tt
        store_in_folder = CmdHistory._storageFolder()
        folder = os.path.join(store_in_folder, tt)
        additional_info_fits_file_name = os.path.join(folder, 'info.fits')
        header = pyfits.getheader(additional_info_fits_file_name)
        hduList = pyfits.open(additional_info_fits_file_name)
        theObject._modeVector = hduList[0].data
        theObject._indexingList = hduList[1].data
        theObject._cmdMatrix = hduList[2].data
        theObject._cmdHToApply = hduList[3].data
        theObject._ampVect = hduList[4].data
        theObject._who = header['WHO']
        try:
            theObject._nPushPull = header['NPUSHPUL']
        except KeyError:
            theObject._nPushPull = 1
        return theObject

    @staticmethod
    def loadFromH5(device, tt):
        """ Creates the object from the info.fits file located in tt

            Args:
                tt = tracking number

            Returns:
                    theObject = CmdHistory class object
        """
        theObject = CmdHistory(device)
        theObject._tt = tt
        store_in_folder = CmdHistory._storageFolder()
        folder = os.path.join(store_in_folder, tt)
        file_name = os.path.join(folder, 'info.h5')
        hf = h5py.File(file_name, 'r')
        hf.keys()
        data1 = hf.get('dataset_1')
        data2 = hf.get('dataset_2')
        data3 = hf.get('dataset_3')
        data4 = hf.get('dataset_4')
        theObject._nPushPull = hf.attrs['NPUSHPUL']
        theObject._who = hf.attrs['WHO']
        theObject._modeVector = np.array(data1)
        theObject._indexingList = np.array(data2)
        theObject._cmdHToApply = np.array(data3)
        theObject._ampVect = np.array(data4)
        hf.close()
        return theObject
