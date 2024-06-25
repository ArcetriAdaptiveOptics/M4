"""
Author(s):
    - P. Ferraiuolo

Written in June 2024
"""
import os
import numpy as np
from m4.configuration import read_iffconfig
from m4.ground import read_data as rd
from m4.configuration import config_folder_names as fn

iffold = fn.IFFUNCTIONS_ROOT_FOLDER

class IFFCapturePreparation():
    """
    Class for the preparation for the Influence Function Acquisition

    Methods
    -------
    _getCmdMatrix

    _createTriggerPadding

    _createRegistrationPattern

    createAuxCmdHistory

    createCmdMatrixHistory

    createTimedCmdHistory
    """
    def __init__(self, dm):
        '''The Constructor'''
        # DM information
        self.mirrorModes        = dm.mirrorModes
        self._NActs             = dm.nActs
        # IFF info
        self.modalBaseId        = None
        self._modesList         = None
        self._modalBase         = self.mirrorModes 
        self._cmdMatrix         = None
        self._indexingList      = None
        self._modesAmp          = None
        self._template          = None
        # Matrices
        self.timedCmdHistory    = None
        self.cmdMatHistory      = None
        self.auxCmdHistory      = None
        self.triggPadCmdHist    = None
        self.regPadCmdHist      = None

    # def initDM_test(self): #used only for debug
    #     dmFold = os.path.join(fn.OPT_DATA_FOLDER,'test')
    #     cmdMatFile =  os.path.join(dmFold,'ff_v_matrix.fits')
    #     hdu = pyfits.open(cmdMatFile)
    #     cmdMat = hdu[0].data
    #     nActs = np.shape(cmdMat)[0]
    #     return cmdMat, nActs

    def createTimedCmdHistory(self,  modesList=None, modesAmp=None, template=None, shuffle=False): 
        """
        Function that creates the final timed command history to be applied

        Parameters
        ----------
        modesList : int | ArrayLike
            List of selected modes to use. Default is None, that means all modes of the base command matrix are used.

        modesAmp : float
            Amplitude of the modes. Default is None, that means the value is loaded from the 'iffconfig.ini' file

        template : int | ArrayLike
            Template for the push-pull measures. List of 1 and -1. Default is None, which means the template is loaded from the 'iffcongig.ini' file.

        shuffle : boolean
            Decide wether to shuffle or not the modes order. Default is False

        Returns
        -------
        timedCmdHist : float | ArrayLike
            Final timed command history, including the trigger padding, the registration pattern and the command matrix history.
        """
        self._modesList = modesList
        self.createCmdMatrixHistory(modesList, modesAmp, template, shuffle)
        self.createAuxCmdHistory()
        cmdHistory = np.hstack((self.auxCmdHistory, self.cmdMatHistory))
        timing = read_iffconfig.getTiming()
        timedCmdHist = np.repeat(cmdHistory, timing, axis=1) 
        self.timedCmdHistory = timedCmdHist
        return timedCmdHist
    
    def getInfoToSave(self):
        """
        Return the data to save as fits files, arranged in a dictionary

        Returns
        -------
        info : dict
            Dictionary containing all the vectors and matrices needed

        """
        info = {'cmdMatrix': self._cmdMatrix,
                'modesList': self._modesList,
                'ampVector': self._modesAmp,
                'indexList': self._indexingList,
                'template' : self._template,
                'shuffle'  : 'To be implemented'
            }
        return info
        

    def createCmdMatrixHistory(self, mlist, modesAmp=None, template=None, shuffle=False):
        """
        Creates the command matrix history for the IFF acquisition.

        Parameters
        ----------
        modesAmp : float
            Amplitude of the modes to be commanded. If no argument is passed, it will be loaded from the configuration file iffConfig.ini
        template : int | ArrayLike
            Template for the push-pull application of the modes. If no argument is passed, it will be loaded from the configuration file iffConfig.ini
        shuffle : boolean
            Decides to wether shuffle or not the order in which the modes are applied. Default is False

        Returns
        -------
        cmd_matrixHistory : float | ArrayLike
            Command matrix history to be applied, with the correct push-pull application, following the desired template.
        """
        if template is None:
            _,_,_, template,_ = read_iffconfig.getConfig('IFFUNC')
        if modesAmp is None:
            _,_,modesAmp,_,_ = read_iffconfig.getConfig('IFFUNC')
        self._template = template
        nModes = self._cmdMatrix.shape[1]
        n_push_pull = len(template)
        if np.size(modesAmp)==1:
            modesAmp = np.full(nModes, modesAmp)
        self._createCmdMatrix(mlist)
        self._modesList = mlist
        self._modesAmp  = modesAmp
        if shuffle:
            cmd_matrix = np.zeros((self._cmdMatrix.shape[0], self._cmdMatrix.shape[1]))
            modesList = np.copy(self._modesList)
            np.random.shuffle(modesList)
            k=0
            for i in modesList:
                cmd_matrix.T[k] = self._cmdMatrix[i]
                k += 1
            self._indexingList = np.array(modesList)
        else:
            cmd_matrix = self._cmdMatrix
            modesList = self._modesList
        n_frame = len(self._modesList) * n_push_pull
        cmd_matrixHistory = np.zeros((self._NActs, n_frame))
        for j in range(n_push_pull):
            for i in range(nModes):
                k = cmd_matrix.shape[1]*j + i
                cmd_matrixHistory.T[k] = cmd_matrix[:,i]*template[j]*modesAmp[i]
        self.cmdMatHistory = cmd_matrixHistory
        return cmd_matrixHistory

    def createAuxCmdHistory(self):
        '''
        Creates the initial part of the final command history matrix that will\
        be passed to M4. This includes the Trigger Frame, the first frame to have\
        a non-zero command, and the Padding Frame, two frames with high rms, useful\
        for setting a start to the real acquisition.

        Result
        ------
        aus_cmdHistory : float | ArrayLike

        '''
        self._createRegistrationPattern()
        self._createTriggerPadding()
        aux_cmdHistory = np.hstack((self.triggPadCmdHist, self.regPadCmdHist))
        self.auxCmdHistory = aux_cmdHistory
        return aux_cmdHistory

    def _createRegistrationPattern(self):
        """
        Creates the registration pattern to apply after the triggering and before\
        the commands to apply for the IFF acquisition. The information about number\
        of zeros, mode(s) and amplitude are read from the 'iffconfig.ini' file.

        Returns
        -------
        regHist : float | ArrayLike
            Registration pattern command history

        """
        nZeros, regId, regAmp, regTemp,mbase = read_iffconfig.getConfig('REGISTRATION')
        self._updateModalBase(mbase)    
        zeroScheme = np.zeros((self._NActs, nZeros))
        regScheme = np.zeros((self._NActs, len(regTemp)*len(regId)))

        k=0
        for i in regId:
            for t in range(len(regTemp)):
                regScheme[i,k] = regAmp*regTemp[t]
                k+=1

        regHist = np.hstack((zeroScheme, regScheme))
        self.regPadCmdHist = regHist

        return regHist

    def _createTriggerPadding(self):
        """
        Function that creates the trigger padding scheme to apply before the \
        registration padding scheme. The information about number of zeros, \
        mode(s) and amplitude are read from the 'iffconfig.ini' file.

        Returns
        -------
        triggHist : float | ArrayLike
            Trigger padding command history

        """
        nZeros, trigId, trigAmp, _, mbase= read_iffconfig.getConfig('TRIGGER')
        self._updateModalBase(mbase)
        zeroScheme = np.zeros((self._NActs, nZeros))
        trigMode = self._modalBase[:,trigId]*trigAmp
        triggHist = np.hstack((zeroScheme, trigMode))
        self.triggPadCmdHist = triggHist
        return triggHist

    def _createCmdMatrix(self, mlist):
        '''
        Cuts the modal base according the given modes list
        '''
        if self.modalBaseId is None:
            _,_,_,_,baseId = read_iffconfig.getConfig('IFFUNC')
        else:
            baseId = self.modalBaseId
        self._updateModalBase(baseId)
        self._cmdMatrix = self._modalBase[:,mlist] 
        return self._cmdMatrix

    def _updateModalBase(self, mbasename: str = None):
        '''
        Redefines the command matrix to be used ma dice marco di buttare
        '''
        print(mbasename)
        if (mbasename is None) or (mbasename == 'mirror'):
            print('Using mirror modes')
            self._modalBase = self.mirrorModes
        elif mbasename == 'zonal':
            print('Using zonal modes')
            self._modalBase = self._createZonalMat()
        elif mbasename == 'hadamard':
            print('Using Hadamard modes')
            self._modalBase = self._createHadamardMat()
            #implement here other options, or they can be read by file as below
        else:
            print('Using user-defined modes')
            self._modalBase = self._createUserMat(mbasename) #this is expected to be a tracknum
            

    def _createUserMat(self, tracknum: str = None):
        """
        

        Parameters
        ----------
        tracknum : str, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        cmdBase : TYPE
            DESCRIPTION.

        """
        print('Reading modal base from tracknum: '+tracknum)
        modalBaseFileName = 'Standard modal base file name' # !!! RE-DEFINE THIS 
        mbfile = os.path.join(fn.MODALBASE_ROOT_FOLDER, tracknum, modalBaseFileName)
        cmdBase = rd.readFits_data(mbfile)
        return cmdBase

    def _createZonalMat(self):
        """
        

        Returns
        -------
        cmdBase : TYPE
            DESCRIPTION.

        """
        cmdBase = np.eye(self._NActs)
        return cmdBase

    def _createHadamardMat(self):
        """
        

        Returns
        -------
        cmdBase : TYPE
            DESCRIPTION.

        """
        from scipy.linalg import hadamard
        import math
        numb = math.ceil(math.log(self._NActs,2))
        hadm = hadamard(2**numb)  #here we suppose that the HADAMARD matrix is used only to measure SEGMENT IFF (which makes sense)
        cmdBase = hadm[1:self._NActs+1, 1:self._NActs+1] #to remove the piston mode. check if it is ok !!!!!
        print('Removed 1st column of Hadamard matrix, or piston mode')
        return cmdBase