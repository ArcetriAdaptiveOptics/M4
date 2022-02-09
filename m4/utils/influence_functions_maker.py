'''
Authors
  - C. Selmi:  written in 2019
'''

import os
import copy
import logging
import h5py
import numpy as np
from astropy.io import fits as pyfits
from m4.ground import tracking_number_folder
from m4.configuration import config_folder_names as fold_name
from m4.type.modalAmplitude import ModalAmplitude
from m4.type.modalBase import ModalBase
from m4.type.commandHistory import CmdHistory


#ad ora pensata per comandare il segmento; manca l'interferometro
class IFFunctionsMaker():
    '''
    This class is responsible for the application of zonal or global commands
    to the deformable mirror in order to collect the functions of influence.
    Data are saving in the folder corresponding to
    the tracking number generated

    HOW TO USE IT::

        creazione oggetto dm e interf
        from m4.influence_functions_maker import IFFunctionsMaker
        IFF = IFFunctionsMaker(dm, interf)
        tt = IFF.acq_IFFunctions(modesVectorFitsFileName, nPushPull,
                                amplitudeFitsFileName, cmdMatrixFitsFileName,
                                shuffle=None, template=None)
    '''

    def __init__(self, deformable_mirror, interf):
        """The constructor """
        self._logger = logging.getLogger('IFF_MAKER:')
        self._nActs = deformable_mirror.getNActs() #OttParameters.N_ACT_SEG
        self._dm = deformable_mirror
        self._interf = interf

        self._nPushPull = None
        self._template = None
        self._amplitudeTag = None
        self._amplitude = None
        self._cmdMatrixTag = None
        self._cmdMatrix = None
        self._actsVector = None
        self._tt_cmdH = None
        self._indexingList = None
        self._tt = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return fold_name.IFFUNCTIONS_ROOT_FOLDER


    def acq_IFFunctions(self, n_push_pull, amplitude_fits_file_name,
                        cmd_matrix_fits_file_name,
                        shuffle=None, template=None):
        '''
        Performs the process of acquiring interferograms

        Parameters
        ----------
             n_push_pull: int
                          number of push pull on the actuator
             amplitude_fits_file_name: string
                                       fits file name
                                       (Example = amp.fits)
             cmd_matrix_fits_file_name: int
                                        fits file name
                                        (Example = modalBase.fits)
        Other Parameters
        ----------------
             shuffle: optional
                      if not indicated, the function create the tidy command
                      history matrix
             template: numpy array, optional
                       vector composed by 1 and -1
                       if not indicated, the function create the vector [1, -1, 1]

        Returns
        -------
                tt: string
                    tracking number of measurements made
        '''
        amplitude, cmd_matrix = self._readTypeFromFitsNameTag(amplitude_fits_file_name,
                                                              cmd_matrix_fits_file_name)

        self._nPushPull = n_push_pull
        if template is None:
            self._template = np.array([1, -1, 1])
        else:
            self._template = template

        self._amplitudeTag = amplitude_fits_file_name
        self._amplitude = amplitude
        self._cmdMatrixTag = cmd_matrix_fits_file_name
        self._cmdMatrix = cmd_matrix

        self._actsVector = np.arange(self._nActs)
        indexing_input = copy.copy(self._actsVector)
        dove, tt = tracking_number_folder.createFolderToStoreMeasurements(self._storageFolder())
        self._tt = tt

#         diagonal_mat = self._diagonalControll(cmd_matrix)
#         if np.count_nonzero(diagonal_mat -
#                             np.diag(np.diagonal(diagonal_mat))) == 0:
#             print('Measure of zonal IF')
#             self._logger.info("Measurement of zonal influence functions for segment?. Location: %s",
#                               tt)
#         else:
#             print('Measure of global IF')
#             self._logger.info("Measurement of modal influence functions for segment?. Location: %s",
#                               tt)


        cmdH = CmdHistory(self._nActs)

        if shuffle is None:
            command_history_matrix_to_apply, self._tt_cmdH = \
                    cmdH.tidyCommandHistoryMaker(self._actsVector,
                                                 amplitude,
                                                 cmd_matrix,
                                                 n_push_pull,
                                                 template)
        else:
            command_history_matrix_to_apply, self._tt_cmdH = \
                    cmdH.shuffleCommandHistoryMaker(self._actsVector,
                                                    amplitude,
                                                    cmd_matrix,
                                                    n_push_pull,
                                                    template)
        self._indexingList = cmdH.getIndexingList()
        #self._saveInfo(dove, 0)
        #self._saveInfoSeparately(dove)

        n_images = 1
        for i in range(command_history_matrix_to_apply.shape[1]):
            for j in range(self._template.shape[0]):
                k = self._template.shape[0]*i +j
                print('Acquisition of image %d' %k)
                self._dm.setActsCommand(command_history_matrix_to_apply[:, i]*j) #vecchia setShape
                masked_image = self._interf.acquire_phasemap(n_images)
                file_name = 'image_%04d.fits' %k
                self._interf.save_phasemap(dove, file_name, masked_image)
        return tt

    def _readTypeFromFitsNameTag(self, amplitude_fits_file_name,
                                 cmd_matrix_fits_file_name):
        '''
        Parameters
        ----------
            amplitude_fits_file_name: string
                                     vector with mode amplitude fits file name
            cmd_matrix_fits_file_name: string
                                    matrix of mode commands fits file name
        Returns
        -------
            amplitude: numpy array
                    vector with mode amplitude
            cmd_matrix: numpy array [nActs x nModes]
                        matrix of mode commands
                        diagonal matrix in case of zonal commands
        '''
        ma = ModalAmplitude.loadFromFits(amplitude_fits_file_name)
        amplitude = ma.getModalAmplitude()

        mb = ModalBase.loadFromFits(cmd_matrix_fits_file_name)
        cmd_matrix = mb.getModalBase()
        return amplitude, cmd_matrix

    def _diagonalControll(self, matrix):
        """ Prepares the matrix for the diagonal control """
        v = np.zeros((self._nActs, 1))
        reps = matrix.shape[0] - matrix.shape[1]
        vects = np.tile(v, reps)
        new_matrix = np.hstack((matrix, vects))
        return new_matrix


    def _testIFFunctions_createCube25fromFileFitsMeasure(self):
        """ Test function: create a measurement cube using data
        generate whit m4 idl software"""
        fold_my_pc = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/OIM_25modes.fits'
        #fold_m4_pc = os.path.join(Configuration.OPD_DATA_FOLDER, 'TestData', 'OIM_25modes.fits')
        hduList = pyfits.open(fold_my_pc)
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


    def _saveInfo(self, folder, fits_or_h5):
        """ Save the fits info file containing the input data for
        the creation of iff

        Args:
            folder = path that indicates where to save the info file
        """
        if fits_or_h5 == 0:
            fits_file_name = os.path.join(folder, 'info.fits')
            header = pyfits.Header()
            header['NPUSHPUL'] = self._nPushPull
            header['TT_CMDH'] = self._tt_cmdH
            header['MODEVECT'] = self._modesVectorTag
            header['CMDMAT'] = self._cmdMatrixTag
            header['AMP'] = self._amplitudeTag
            pyfits.writeto(fits_file_name, self._indexingList, header)
            pyfits.append(fits_file_name, self._template)
        else:
            fits_file_name = os.path.join(folder, 'info.h5')
            hf = h5py.File(fits_file_name, 'w')
            hf.create_dataset('dataset_1', data=self._indexingList)
            hf.attrs['MODEVECT'] = self._modesVectorTag
            hf.attrs['CMDMAT'] = self._cmdMatrixTag
            hf.attrs['AMP'] = self._amplitudeTag
            hf.attrs['NPUSHPUL'] = self._nPushPull
            hf.attrs['TT_CMDH'] = self._tt_cmdH
            hf.close()

    def _saveInfoSeparately(self, folder):
        """ Save the input data for the creation of iff

        Args:
            folder = path that indicates where to save
        """
        fits_file_name = os.path.join(folder, 'indexingList.fits')
        pyfits.writeto(fits_file_name, self._indexingList)
        fits_file_name = os.path.join(folder, 'template.fits')
        pyfits.writeto(fits_file_name, self._template)
        fits_file_name = os.path.join(folder, 'more_info.txt')
        file = open(fits_file_name, 'w+')
        file.write('tt_cmdH = %s' %(self._tt_cmdH))
        file.close()
        fits_file_name = os.path.join(folder, 'tag_info.txt')
        file = open(fits_file_name, 'w+')
        file.write('Modes_vector_tag = %s, Cmd_matrix_tag = %s, Amplitude_tag = %s' \
                   %(self._modesVectorTag, self._cmdMatrixTag, self._amplitudeTag))
        file.close()
        fits_file_name = os.path.join(folder, 'n_push_pull.txt')
        file = open(fits_file_name, 'w+')
        file.write('N_push_pull = %e' %(self._nPushPull))
        file.close()

    @staticmethod
    def loadInfo(tt, fits_or_h5=0):
        """ Reload information contained in fits Info file

        Parameters
        ----------
        fits_file_name: string
                        info file path

        Returns
        -------
        who: string
            which segment
        tt_cmdH: string
                CommandHistory tracking number
        acts_vector: numpy array
                    vector of actuators
        cmd_matrix: matrix
                    modal base
        cmd_ampl: numpy array
                amplitude vector
        n_push_pull: int
                    number of push pull
        indexingList: list
                    list of index order using in command history
        """
        store_in_folder = IFFunctionsMaker._storageFolder()
        folder = os.path.join(store_in_folder, tt)
        dm = 0
        interf = 0
        theObject = IFFunctionsMaker(dm, interf)
        if fits_or_h5 == 0 :
            additional_info_fits_file_name = os.path.join(folder, 'info.fits')
            header = pyfits.getheader(additional_info_fits_file_name)
            hduList = pyfits.open(additional_info_fits_file_name)
            theObject._modesVectorTag = header['MODEVECT']
            theObject._cmdMatrixTag = header['CMDMAT']
            theObject._amplitudeTag = header['AMP']
            theObject._indexingList = hduList[0].data
            theObject._template = hduList[1].data
            theObject._tt_cmdH = header['TT_CMDH']
            try:
                theObject._nPushPull = header['NPUSHPUL']
            except KeyError:
                theObject._nPushPull = 1
        else:
            file_name = os.path.join(folder, 'info.h5')
            hf = h5py.File(file_name, 'r')
            hf.keys()
            data_1 = hf.get('dataset_1')
            theObject._modesVectorTag = hf.attrs['MODEVECT']
            theObject._cmdMatrixTag = hf.attrs['CMDMAT']
            theObject._amplitudeTag = hf.attrs['AMP']
            theObject._indexingList = np.array(data_1)
            theObject._nPushPull = hf.attrs['NPUSHPUL']
            theObject._tt_cmdH = hf.attrs['TT_CMDH']
        theObject._h5Folder = folder
        return theObject
