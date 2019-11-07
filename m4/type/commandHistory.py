'''
@author: cs
'''
import os
import copy
import pyfits
import h5py
import numpy as np
from m4.ground import tracking_number_folder
from m4.ground.configuration import Configuration
from m4.ground import logger



class CmdHistory():
    '''
    Classe per creare e gestire la matrice con la storia dei comandi
    '''

    def __init__(self, device):
        self._device = device
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
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                            "CommandHistory")
    '''
         arg:
             modesVector= vettore dell'indice dei modi o degli attuatori
                        da applicare (numpy.array([]))
             nPushPull= numero di push pull consecutivi sull'attuatore
                         (int)
             ampVector= vettore con l'ampiezza dei modi (numpy.array([]))
             cmdMatrix= matrice dei comandi dei modi 
                         (nActs x nModes)
                         matrice diagonale nel caso di comandi zonali
         return:
             matrixToApply: command history (nAct, nModes x nPushPusll x 2)
    '''

    def tidyCommandHistoryMaker(self, mode_vector, amp_vector,
                                cmd_matrix, n_push_pull):
        self._ampVect = amp_vector
        cmd_history = self._tidyCmdHistory(mode_vector, n_push_pull, cmd_matrix)
        aa = np.arange(self._cmdHistory.shape[1])
        bb = np.tile(amp_vector, n_push_pull)
        zipped = zip(aa, bb)
        matrix_to_apply = self._cmdHistoryToApply(zipped)
        self._cmdHToApply = matrix_to_apply
        tt = self.saveAsFits()
        logger.log('Creation of the ordered', 'commandHistoryMatrix', tt)
        print(tt)

        return matrix_to_apply, tt


    def shuffleCommandHistoryMaker(self, mode_vector, amp_vector,
                                   cmd_matrix, n_push_pull):
        self._ampVect = amp_vector
        cmd_history, indexing_list = self._shuffleCmdHistory(
            mode_vector, n_push_pull, cmd_matrix)
        zipped = self._zippedAmplitude(amp_vector)
        matrix_to_apply = self._cmdHistoryToApply(zipped)
        self._cmdHToApply = matrix_to_apply
        tt = self.saveAsFits()
        logger.log('Creation of the shuffle', 'commandHistoryMatrix', tt)
        print(tt)

        return matrix_to_apply, tt

    def _shuffleCmdHistory(self, mode_vector, n_push_pull, cmd_matrix):
        self._modeVector = copy.copy(mode_vector)
        self._nPushPull = n_push_pull
        self._cmdMatrix = cmd_matrix

        n_frame = mode_vector.size * n_push_pull
        matrix_to_apply = np.zeros((self._nActs, n_frame))

        indexingList = []
        for j in range(n_push_pull):
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


    def _tidyCmdHistory(self, mode_vector, n_push_pull, cmd_matrix):
        self._modeVector = copy.copy(mode_vector)
        self._nPushPull = n_push_pull
        self._cmdMatrix = cmd_matrix
        indList = []
        for i in range(n_push_pull):
            indList.append(mode_vector)
        self._indexingList = np.array(indList)
        #self._indexingList= np.tile(mode_vector, n_push_pull)

        n_frame = mode_vector.size * n_push_pull
        matrix_to_apply = np.zeros((self._nActs, n_frame))

        cmdList = []
        for i in mode_vector:
            cmd = cmd_matrix[:, i]
            cmdList.append(cmd)

        for j in range(n_push_pull):
            for i in range(len(cmdList)):
                k = len(cmdList)*j + i
                matrix_to_apply.T[k] = cmdList[i]

        self._cmdHistory = matrix_to_apply

        return self._cmdHistory


    def _amplitudeReorganization(self, indexing_input, indexing_list,
                                 amplitude, n_push_pull):
        where = []
        for i in indexing_input:
            for j in range(n_push_pull):
                a = np.where(indexing_list[j] == i)
                where.append(a)


        where = np.array(where)
        vect = np.zeros(amplitude.shape[0]*n_push_pull)

        for i in range(amplitude.shape[0]):
            for k in range(n_push_pull):
                p = n_push_pull * i + k
                indvect = where[p][0][0]+ indexing_input.shape[0] * k
                vect[indvect] = amplitude[i]

        return vect


    def _zippedAmplitude(self, amplitude):
        aa = np.arange(self._cmdHistory.shape[1])
        reorganized_amplitude = \
                self._amplitudeReorganization(self._modeVector,
                                              self._indexingList,
                                              amplitude,
                                              self._nPushPull)
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
