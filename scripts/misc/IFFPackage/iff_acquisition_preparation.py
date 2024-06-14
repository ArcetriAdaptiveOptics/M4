'''
NOTA:
    la classe ha la creazione/lettura della base modale e i metodi per farci le storie temporali
    la classe è replicata, ereditando i metodi, e aggiornando la crezione della base modale secondo cosa richiesto
    esempio: istanzio la classe con la base modale che voglio e la classe si è ereditata i metodi della classe base
'''
'''
Author(s):
    - P. Ferraiuolo

Written in June 2024
'''
'''
WARNING!!! al lthe elements relative to a cmdHistory (trigger, padding, etc), shall be act x frames (e.g. 892 x 200 commands)

'''
import numpy as np
import os
from astropy.io import fits as pyfits
from m4.mini_OTT import timehistory as th
from m4.configuration import read_iffconfig
from m4.configuration import config_folder_names as fn
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
    def __init__(self): #file .ini per inizializzare?
        '''The Constructor'''
        #self._NActs             = read_iffconfig.getNActs_fromConf() #Dove trovo l'info? Si deve caricare il dm? o leggiamo conf?
        self._modesList         = None
        #self._cmdMatrix=None #        = initDM_test()
        self._mirrorModalBase, self._NActs = self.initDM_test()
        self._modalBase         = self._mirrorModalBase   #e poi modificare questo quando si ridefinisce la base modale (zonal, had, ...)
        self._cmdMatrix         = None

        self._indexingList      = None

        self.timedCmdHistory    = None

        self.cmdMatHistory      = None
        self.auxCmdHistory      = None

        self.triggPadCmdHist    = None
        self.regPadCmdHist      = None

    def initDM_test(self):
        '''
        This function serves only for debugging, to retrieve the DM data
        '''
        dmFold = os.path.join(fn.OPT_DATA_FOLDER,'test')
        cmdMatFile =  os.path.join(dmFold,'ff_v_matrix.fits')
        hdu = pyfits.open(cmdMatFile)
        cmdMat = hdu[0].data
        nActs = np.shape(cmdMat)[0]
        return cmdMat, nActs

    def createTimedCmdHistory(self, cmdBase, modesList=None, modesAmp=None, template=None, shuffle=False):
        """
        Function that creates the final timed command history to be applied

        Parameters
        ----------
        cmdBase : str
           Identification string of the base nd matrix to load. Values can be:

                'zonal'    : It loads the zonal command matrix of the mirror;

                'hadamard' : It loads the hadamard's command matrix;

                'tn'       : A tracking number in which is contained the fits file of the command matrix to load. Can be the case of the mirror's modal command matrix.

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
        self._getCmdMatrix(cmdBase, modesList)
        self.createCmdMatrixHistory(modesAmp, template, shuffle)
        self.createAuxCmdHistory()
        cmdHistory = np.hstack((self.auxCmdHistory, self.cmdMatHistory))

        timing = read_iffconfig.getTiming()

        timedCmdHist = np.repeat(cmdHistory, timing, axis=1) # Timing info where? iffconfig.ini?

        self.timedCmdHistory = timedCmdHist
        return timedCmdHist

    def createCmdMatrixHistory(self, modesAmp=None, template=None, shuffle=False):
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
            _,_,_, template = read_iffconfig.getConfig('IFFUNC')
        if modesAmp is None:
            _,_,modesAmp,_ = read_iffconfig.getConfig('IFFUNC')

        n_push_pull = len(template)

        if shuffle:
            cmd_matrix = np.zeros((self._cmdMatrix.shape[0], self._cmdMatrix.shape[1]))
            modesList = np.random.shuffle(self._modesList)
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
            for i in range(cmd_matrix.shape[1]):
                k = cmd_matrix.shape[1]*j + i
                cmd_matrixHistory.T[k] = cmd_matrix[:,i]*template[j]*modesAmp

        self.cmdMatHistory = cmd_matrixHistory
        return cmd_matrixHistory

    def createAuxCmdHistory(self):
        '''
        Creates the initial parteof the final command history matrix that will be passed to M4. This includes the Trigger Frame, the first frame to have a non-zero command, and the Padding Frame, two frames with high rms, useful for setting a start to the real acquisition.

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
        Creates the registration pattern to apply after the triggering and before the commands to apply for the IFF acquisition. The information about number of zeros, mode(s) and amplitude are read from the 'iffconfig.ini' file.

        Returns
        -------
        regHist : float | ArrayLike
            Registration pattern command history

        """
        nZeros, regId, regAmp, regTemp = read_iffconfig.getConfig('REGISTRATION')

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
        Function that creates the trigger padding scheme to apply before the registration padding scheme. The information about number of zeros, mode(s) and amplitude are read from the 'iffconfig.ini' file.

        Returns
        -------
        triggHist : float | ArrayLike
            Trigger padding command history

        """
        nZeros, trigId, trigAmp,_= read_iffconfig.getConfig('TRIGGER')
        zeroScheme = np.zeros((self._NActs, nZeros))
        trigMode = self._mirrorModalBase[:,trigId]*trigAmp    #modRB

        #modRB
        #if self._cmdMatrix is not None:
        #    trigMode = self._modalBase[:,trigId]*trigAmp
        #else: raise ValueError("Command Matrix not Loaded yet. Run 'createCmdMatrixHistory' first to solve the issue")

        triggHist = np.hstack((zeroScheme, trigMode))
        self.triggPadCmdHist = triggHist

        return triggHist

    def _createCmdMatrix(self,mlist):
        '''
        Cuts the modal base according the given modes list
        '''
        self._cmdMatrix = self._modalBase[:,mlist] #qui chiarire la direzione. è così per i dati in arrivo da IDL


    def _getCmdMatrix2(self, identif: str, mlist: list = None):  #originale di Pietro
        """
        This function gets the base command matrix to use for the IFF acquisition. If 'zonal' or 'hadamard' are passe as arguments, it is analytically generated, while if a tn is passed, it will load the modal base from the corresponding .fits file contained in the tn folder.

        Parameters
        ----------
        identif : str
            Identification string of the command matrix to load. Values can be:

                'zonal'    : It loads the zonal command matrix of the mirror;

                'hadamard' : It loads the hadamard's command matrix;

                'tn'       : A tracking number in which is contained the fits file of the command matrix to load. Can be the case of the mirror's modal command matrix.

        mlist : int | int or ArrayLike
            List of modes to be used from the loaded command matrix.

        Returns
        -------
        cmdMat : float | ArrayLike
            Extracted command matrix of the selected modes only.

        """
        if isinstance(identif, str) is True:
            if identif=='zonal':
                cmdBase = np.eye(self._NActs)
            elif identif=='hadamard':
                from scipy.linalg import hadamard
                hadm = hadamard(2**10)  #here we suppose that the HADAMARD matrix is used only to measure SEGMENT IFF (which makes sense)
                cmdBase = hadm[:self._NActs, :self._NActs]
            else:
                flist = th.fileList(identif)
                for item in flist:
                    if 'nomefile' in item:
                        file = item
                with pyfits.open(file) as hdul:
                    cmdBase = hdul[0].data
        else:
            raise TypeError("'identif' must be a str. Accepted values are 'zonal', 'hadamard' or a tracking number")

        if mlist is not None:
            self._indexingList = np.linspace(0, len(mlist), 1)
            cmdMat = np.zeros((cmdBase.shape[0], len(mlist)))
            k = 0
            for n in mlist:
                cmdMat.T[k] = cmdBase[:,n]
                k += 1

            self._cmdMatrix = cmdMat
        else:
            self._cmdMatrix = cmdBase
            self._indexingList = self._modesList = np.arange(0, cmdBase.shape[1], 1)

        return self._cmdMatrix

    def updateCommandMatrix(self, identif: str =None):
        '''
        Redefines the command matrix to be used ma dice marco di buttare
        '''
        if isinstance(identif, str) is True:
            if identif=='zonal':
                cmdBase = np.eye(self._NActs)
            elif identif=='hadamard':
                from scipy.linalg import hadamard
                hadm = hadamard(2**10)  #here we suppose that the HADAMARD matrix is used only to measure SEGMENT IFF (which makes sense)
                cmdBase = hadm[:self._NActs, :self._NActs]
            else:
                flist = th.fileList(identif)
                for item in flist:
                    if 'nomefile' in item:
                        file = item
                with pyfits.open(file) as hdul:
                    cmdBase = hdul[0].data
        else:
            cmdBase = self._cmdMatrix
            print('Using CmdMatrix as passed from mirror configuration, nothing done')
            #raise TypeError("'identif' must be a str. Accepted values are 'zonal', 'hadamard' or a tracking number")


    def _saveMatrix(filename,matrix):
        '''
        Function as a placeholder to remember to save the matrices
        !!! warning! shall this be inside a class? or shall it be a static method?
        we don't want to create the object (== updating the intmat or doing anything weird) just to save the file
        '''
        pass
