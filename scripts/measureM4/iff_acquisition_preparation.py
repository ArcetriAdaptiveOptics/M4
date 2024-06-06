'''
Author(s):
    - P. Ferraiuolo

Written in June 2024
'''
import numpy as np
from astropy.io import fits as pyfits
from m4.mini_ott import timehisory as th
from m4.configuration import iffConfig as iffc  #contains a dictionary or similar with configuration

class IFFCapturePreparation():

    def __init__(self): #file .ini per inizializzare? No, non c'è bisogno. leggono le funzioni in caso
        self._nActs             = ifcc.NACTS
        self._cmdMatrix         = None
        self._indexingList      = None
        self._modeVector        = None
        
        self.cmdMatHistory      = None
        self.auxCmdHistory      = None
        self.triggPadCmdList    = None
        self.regPadCmdList      = None


    def createTimedCmdHistory(self): 
        '''
        Function that creates the final timed command history to be applied

        Parameters
        ----------

        Returns
        -------
        timedCmdHist : float | ArrayLike
            Final timed command history, including the trigger padding, the registration pattern and the command matrix history.
        '''
        self.createCmdMatrixHistory()
        self.creatAuxCmdHistory()
        cmdHistory = np.hstack((slef.auxCmdHistory, self.cmdMatHistory))

        timedCmdHist = np.repeat((cmdHistory, ifcc.TIMING, axis=1)

        self.triggPadCmdHistory = timedCmdhist
        return timedCmdHist

    def createCmdMatrixHistory(self, modesList, modes_amplitude, template=None, shuffle=False):
        '''

        '''
        if template is None:
            template = [1,-1,1]

        self._modeVector = modesList
        n_push_pull      = len(template)
        cmd_matrix       = self._getCmdMatrix(ifcc.CMD_BASE, modesList)

        indList = []
        if shuffle==True:
            for i in range(len(n_push_pull):
                np.random.shuffle(modesList)
                indexingList.append(list(modesList))
        else:
            for i in range(n_push_pull):
                indList.append(modesList)

        self._indexingList = np.array(indList)

        n_frame = len(modesList) * n_push_pull
        cmd_matrixHistroy = np.zeros((self._nActs, n_frame))

        cmdList = []
        for i in mode_vector:
            cmd = cmd_matrix[:, i]
            cmdList.append(cmd)

        for j in range(n_push_pull):
            for i in range(len(cmdList)):
                k = len(cmdList)*j + i
                cmd_matrixHistory.T[k] = cmdList[i]

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

    def _createRegistrationPattern(self, reg_modes, reg_amplitudes, reg_template):
        '''
        Creates the registration pattern to apply after the triggering and before the commands to apply for the IFF acquisition

        Returns
        -------
        regHist : float | ArrayLike
            Registration pattern command history

        '''

        #stesso schema di trigger padding
        #input viene da iffc: e sarebbe quali attuatori alzare per fare i marker, con che ampiezza ed eventualmente quante ripetizioni (template)

        return regHist

    def _createTriggerPadding(self):
        '''
        Function that creates the trigger padding scheme to apply before the registration padding scheme

        Returns
        -------
        triggHist : float | ArrayLike
            Trigger padding command history

        '''
        paddingScheme = np.zeros((self._NActs, ifcc.N_PADDING_ZEROS))
        trigg = iffc.T0  # trigg o è un vettore con già le ampiezze desiderate, oppure si crea un vettore vuoto e si popola. Meglio la prima se è roba di default sempre uguale
        triggHist = np.hstack((paddingScheme, trigg))
        
        self.triggPadCmdList = triggHist

        return triggHist


    def _getCmdMatrix(identif: str, mlist: list):
        """
        This function gets the base command matrix to use for the IFF acquisition. If 'zonal' or 'hadamard' are passe as arguments, it is analytically generated, while if a tn is passed, it will load the modal base from the corresponding .fits file contained in the tn folder.

        Parameters
        ----------
        identif : str
            Identification string of the command matrix to load. Values can be:

                'zonal'    : It loads the zonal command matrix of the mirror;

                'hadamard' : It loads the hadamard's command matrix;

                'tn'       : A tracking number containing the fits file of the command matrix to load. Can be the case of the mirror's modal command matrix.

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

    cmdMat = np.zeros((cmdBase.shape[0], len(mlist)))
    k = 0
    for n in mlist:
        cmdMat.T[k] = cmdBase[:,n]
        k += 1

    self._cmdMatrix = cmdMat
    return cmdMat
