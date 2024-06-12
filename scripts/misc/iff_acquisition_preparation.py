'''
Author(s):
    - P. Ferraiuolo

Written in June 2024
'''
import numpy as np
from astropy.io import fits as pyfits
from m4.mini_ott import timehisory as th
from m4.configuration.read_iffconfig import getConfig

bpath = os.environ['M4CONF']

class IFFCapturePreparation():

    def __init__(self): #file .ini per inizializzare?
        '''The Constructor'''
        self._nActs             = ifcc.NACTS #Dove trovo l'info? Si deve caricare il dm? o leggiamo conf?
        self._modesList         = None
        self._cmdMatrix         = None
        self._indexingList      = None
        
        self.cmdMatHistory      = None
        self.auxCmdHistory      = None
        self.triggPadCmdHist    = None
        self.regPadCmdHist      = None


    def createTimedCmdHistory(self, cmdBase, modesList=None, modesAmp=None, template=None, shuffle=False): 
        '''
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
        '''
        self._modesList = modesList
        self._getCmdMatrix(cmdBase, modesList)
        self.createCmdMatrixHistory(modesAmp, template, shuffle)
        self.creatAuxCmdHistory()
        cmdHistory = np.hstack((slef.auxCmdHistory, self.cmdMatHistory))

        timedCmdHist = np.repeat((cmdHistory, ifcc.TIMING, axis=1)) # Timing info where? iffconfig.ini?

        self.triggPadCmdHistory = timedCmdhist
        return timedCmdHist

    def createCmdMatrixHistory(self, modesAmp=None, template=None, shuffle=False):
        '''
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
        '''
        if template is None:
            _,_,_, template = getconfig(bpath, 'IFFUNC')
        if modesAmps is None:
            _,_,modesAmp,_ = getconfig(bpath, 'IFFUNC')

        n_push_pull      = len(template)

        if shuffle:
            cmd_matrix = np.zeros((self._cmdMatrix.shape[0], self._cmdMatrix.shape[1]))
            modesList = np.random.shuffle(self._modesList)
            k=0
            for i in randomList:
                cmd_matrix.T[k] = self._cmdMatrix[i]
                k += 1
            self._indexingList = np.array(modesList)
        else: 
            cmd_matrix = self._cmdMatrix
            modesList = self._modesList

        n_frame = len(self._modesList) * n_push_pull
        cmd_matrixHistroy = np.zeros((self._nActs, n_frame))

        cmdList = []
        for i in modesList:
            cmd = cmd_matrix[:, i]
            cmdList.append(cmd)

        for j in range(n_push_pull):
            for i in range(len(cmdList)):
                k = len(cmdList)*j + i
                cmd_matrixHistory.T[k] = cmdList[i]*template[j]*modesAmp

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

        aux_cmdHistory = np.hstack((self._triggPadCmdHist, self._regPadCmdHist))
        self._auxCmdHistory = aux_cmdHistory

        return aux_cmdHistory

    # Ha senso tenere le due funzioni sotto separate da questa sopra? Potrebbee non essere necessario

    def _createRegistrationPattern(self):
        """
        Creates the registration pattern to apply after the triggering and before the commands to apply for the IFF acquisition. The information about number of zeros, mode(s) and amplitude are read from the 'iffconfig.ini' file.

        Returns
        -------
        regHist : float | ArrayLike
            Registration pattern command history

        """
        nZeros, regId, regAmp, regTemp = getconfig(bpath, 'REGISTRATION')

        zeroScheme = np.zeros((self._NActs, nZeros))
        regScheme = np.zeros((self._NActs, len(regTemp)))

        for t in range(len(regTemp)):
            for i in redIg:
                regScheme[i,t] = regAmp*regTemp[t]

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
        nZeros, trigId, trigAmp, trigTemp = getConfig(bpath, 'TRIGGER')
        zeroScheme = np.zeros((self._NActs, nZeros))

        if self.cmdMatrix is not None:
            trigMode = self.cmdMatrix[:,trigId]*trigAmp
        else raise ValueError("Command Matrix not Loaded yet. Run 'createCmdMatrixHistory' first to solve the issue ")

        triggHist = np.hstack((zeroScheme, trigMode))
        self.triggPadCmdHist = triggHist

        return triggHist


    def _getCmdMatrix(identif: str, mlist: list = None):
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
    if isinstance(indentif, str):
        if identif=='zonal':
            cmdBase = np.eye(self._NActs)
        elif indetif=='hadamard':
            from scipy.linalg import hadamard
            hadm = hadamard(2*10)
            cmdBase = hadm[:self._NActs, :self._NActs]
        else:
            flist = th.fileList(identif)
            for item in flist:
                if 'nomefile' in item:
                    file = item

            with pyfits.open(file) as hdul:
                cmdBase = hdul[0].data
    else:
        raise TypeError("'identif' must be a str, and can be 'zonal', 'hadamard' or a tracking number")
    
    if mlist is not None:
        self._indexingList = np.linspace(0, len(mlist), 1)
        cmdMat = np.zeros((cmdBase.shape[0], len(mlist)))
        k = 0
        for n in mlist:
            cmdMat.T[k] = cmdBase[:,n]
            k += 1

        self._cmdMatrix = cmdMat
        return cmdMat
    else:
        self.cmdMatrix = cmdBase
        self._indexingList = np.linspace(0, cmdBas.shape[1], 1)
        return cmdBase
