"""
Author(s):
----------
    - Pietro Ferraiuolo

Written in June 2024

Description
-----------
This module contains the IFFCapturePreparation class, a class which serves as a
preparator for the Influence Function acquisition by M4, creating the timed com
mand matrix history that will be ultimately used.
More information on its use can be found on the class documentation.
"""
import os
import numpy as np
from m4.configuration import read_iffconfig, config_folder_names as fn
from m4.ground import read_data as rd
from m4.dmutils.iff_processing import _getAcqInfo
iffold = fn.IFFUNCTIONS_ROOT_FOLDER

class IFFCapturePreparation():
    """
    Class containing all the functions necessary to create the final timed
    command matrix history to be executed by M4

    Import and Initialization
    -------------------------
    Import the module and initialize the class with a deformable mirror object

    >>> from m4.dmutils.iff_acquisition_preparation import IFFCapturePreparation
    >>> from m4.devices import deformable_mirror as dm
    >>> m4u = dm.M4AU()
    >>> ifa = IFFCapturePreparation(m4u)

    Methods
    -------
    createTimedCmdHistory

        Creates the final timed command matrix history. Takes 4 positional optional
        arguments, which will be read from a configuration file if not passed

    createCmdMatrixhistory

        Takes the modal base loaded into the class (which can be updated using
        the sub-method _updateModalBase) and returns the wanted command matrix
        with the dedired modes and amplitudes, which can be either passed on as
        arguments or read automatically from a configuration file.

        >>> # As example, wanting to update the modal base using a zonal one
        >>> ifa._updateModalBase('zonal')
        'Using zonal modes'

    createAuxCmdHistory

        Creates the auxiliary command matrix to attach to the command matrix
        history. This auxiliary matrix comprehends the trigger padding and the
        registration padding schemes. the parameters on how to create these
        schemes is written in a configuration file.

    getInfoToSave

        A function that returns a dictionary containing all the useful information
        to save, such as the command matrix used, the used mode list, the indexing
        the amplitudes, the used tamplate and the shuffle option.

    Notes
    -----
    In order for the module to work properly, the tower initialization must be
    run, so that the folder names configuration file is populated.
    From the IPython console

    >>> run '/path/to/m4/initOTT.py'
    >>> from m4.dmutils import iff_acquisition_preparation

    At this point you can either use the dm instance already present in the ran
    file, most likely making the IFFCapturePreparation class to use a FakeDM to
    initialize (might not work), or define a second dm instance

    >>> from m4.devices import deformable_mirror as dfm
    >>> ifa = iff_acquisition_preparation.IFFCapturePreparation(dfm.M4AU())

    Upon developing the deformable_mirror module, the initialization issue will
    be addressed.
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
        self._regActs           = None
        self._cmdMatrix         = None
        self._indexingList      = None
        self._modesAmp          = None
        self._template          = None
        self._shuffle           = 0
        # Matrices
        self.timedCmdHistory    = None
        self.cmdMatHistory      = None
        self.auxCmdHistory      = None
        self.triggPadCmdHist    = None
        self.regPadCmdHist      = None

    def createTimedCmdHistory(self,  modesList=None, modesAmp=None, template=None, shuffle=False):

        """
        Function that creates the final timed command history to be applied

        Parameters
        ----------
        modesList : int | ArrayLike
            List of selected modes to use. Default is None, that means all modes
            of the base command matrix are used.
        modesAmp : float
            Amplitude of the modes. Default is None, that means the value is
            loaded from the 'iffconfig.ini' file
        template : int | ArrayLike
            Template for the push-pull measures. List of 1 and -1. Default is
            None, which means the template is loaded from the 'iffcongig.ini' file.
        shuffle : boolean
            Decide wether to shuffle or not the modes order. Default is False

        Returns
        -------
        timedCmdHist : float | ArrayLike
            Final timed command history, including the trigger padding, the
            registration pattern and the command matrix history.
        """
        self.createCmdMatrixHistory(modesList, modesAmp, template, shuffle)
        self.createAuxCmdHistory()
        if not self.auxCmdHistory is None:
            cmdHistory = np.hstack((self.auxCmdHistory, self.cmdMatHistory))
        else:
            cmdHistory = self.cmdMatHistory
            self._regActs = np.array([])
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
        info = {'cmdMatrix'  : self._cmdMatrix,
                'modesVector': self._modesList,
                'regActs'    : self._regActs,
                'ampVector'  : self._modesAmp,
                'indexList'  : self._indexingList,
                'template'   : self._template,
                'shuffle'    : self._shuffle
            }
        return info

    def createCmdMatrixHistory(self, mlist=None, modesAmp=None, template=None, shuffle=False):

        """
        Creates the command matrix history for the IFF acquisition.

        Parameters
        ----------
        modesAmp : float
            Amplitude of the modes to be commanded. If no argument is passed,
            it will be loaded from the configuration file iffConfig.ini
        template : int | ArrayLike
            Template for the push-pull application of the modes. If no argument
            is passed, it will be loaded from the configuration file iffConfig.ini
        shuffle : boolean
            Decides to wether shuffle or not the order in which the modes are
            applied. Default is False

        Returns
        -------
        cmd_matrixHistory : float | ArrayLike
            Command matrix history to be applied, with the correct push-pull
            application, following the desired template.
        """
        _,_,infoIF = _getAcqInfo()
        if mlist is None:
            mlist = infoIF.get('modes')
        else:
            mlist = mlist
            infoIF['modes'] = mlist
            read_iffconfig.updateConfigFile('IFFUNC', infoIF)
        self._modesList = mlist
        modesAmp = modesAmp if modesAmp is not None else infoIF.get('amplitude')
        template = template if template is not None else infoIF.get('template')
        zeroScheme = infoIF['zeros']
        self._template = template
        self._createCmdMatrix(mlist)
        nModes = self._cmdMatrix.shape[1]
        n_push_pull = len(template)
        if np.size(modesAmp)==1:
            modesAmp = np.full(nModes, modesAmp)
        self._createCmdMatrix(mlist)
        self._modesList = mlist
        self._modesAmp  = modesAmp
        if shuffle is not False:
            self._shuffle = shuffle
            cmd_matrix = np.zeros((self._cmdMatrix.shape[0], self._cmdMatrix.shape[1]))
            modesList = np.copy(self._modesList)
            np.random.shuffle(modesList)
            k=0
            for i in modesList:
                cmd_matrix.T[k] = self._cmdMatrix[i]
                k += 1
            self._indexingList = np.arange(0, len(modesList), 1)
        else:
            cmd_matrix = self._cmdMatrix
            modesList = self._modesList
            self._indexingList = np.arange(0, len(modesList), 1)
        n_frame = len(self._modesList) * n_push_pull
        cmd_matrixHistory = np.zeros((self._NActs, n_frame+zeroScheme))
        k = zeroScheme
        for i in range(nModes):
            for j in range(n_push_pull):
                # k = zeroScheme + cmd_matrix.shape[1]*j + i
                cmd_matrixHistory.T[k] = cmd_matrix[:,i]*template[j]*modesAmp[i]
                k += 1
        self.cmdMatHistory = cmd_matrixHistory
        return cmd_matrixHistory

    def createAuxCmdHistory(self):

        '''
        Creates the initial part of the final command history matrix that will
        be passed to M4. This includes the Trigger Frame, the first frame to
        have a non-zero command, and the Padding Frame, two frames with high
        rms, useful for setting a start to the real acquisition.

        Result
        ------
        aus_cmdHistory : float | ArrayLike

        '''
        self._createTriggerPadding()
        self._createRegistrationPattern()
        if self.triggPadCmdHist is not None and self.regPadCmdHist is not None:
            aux_cmdHistory = np.hstack((self.triggPadCmdHist, self.regPadCmdHist))
        elif self.triggPadCmdHist is not None:
            aux_cmdHistory = self.triggPadCmdHist
        elif self.regPadCmdHist is not None:
            aux_cmdHistory = self.regPadCmdHist
        else:
            aux_cmdHistory = None
        self.auxCmdHistory = aux_cmdHistory
        return aux_cmdHistory

    def _createRegistrationPattern(self):

        """
        Creates the registration pattern to apply after the triggering and before
        the commands to apply for the IFF acquisition. The information about number
        of zeros, mode(s) and amplitude are read from the 'iffconfig.ini' file.

        Returns
        -------
        regHist : float | ArrayLike
            Registration pattern command history

        """
        infoR = read_iffconfig.getConfig('REGISTRATION')
        if len(infoR['modes']) == 0:
            self._regActs = infoR['modes']
            return
        self._regActs = infoR['modes']
        self._updateModalBase(infoR['modalBase'])
        zeroScheme = np.zeros((self._NActs, infoR['zeros']))
        regScheme = np.zeros((self._NActs, len(infoR['template'])*len(infoR['modes'])))
        k=0
        for mode in infoR['modes']:
            for t in range(len(infoR['template'])):
                regScheme.T[k] = self._modalBase.T[mode]*infoR['amplitude']*infoR['template'][t]
                k+=1
        regHist = np.hstack((zeroScheme, regScheme))
        self.regPadCmdHist = regHist
        return regHist

    def _createTriggerPadding(self):

        """
        Function that creates the trigger padding scheme to apply before the
        registration padding scheme. The information about number of zeros,
        mode(s) and amplitude are read from the 'iffconfig.ini' file.

        Returns
        -------
        triggHist : float | ArrayLike
            Trigger padding command history
        """
        infoT = read_iffconfig.getConfig('TRIGGER')
        if len(infoT['modes']) == 0:
            return
        self._updateModalBase(infoT['modalBase'])
        zeroScheme = np.zeros((self._NActs, infoT['zeros']))
        trigMode = self._modalBase[:,infoT['modes']]*infoT['amplitude']
        triggHist = np.hstack((zeroScheme, trigMode))
        self.triggPadCmdHist = triggHist
        return triggHist

    def _createCmdMatrix(self, mlist):

        '''
        Cuts the modal base according the given modes list
        '''
        infoIF = read_iffconfig.getConfig('IFFUNC')
        self._updateModalBase(infoIF['modalBase'])
        self._cmdMatrix = self._modalBase[:,mlist]
        return self._cmdMatrix

    def _updateModalBase(self, mbasename: str = None):

        """
        Updates the used modal base

        Parameters
        ----------
        mbasename : str, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        if (mbasename is None) or (mbasename == 'mirror'):
            #print('Using mirror modes')
            self.modalBaseId = mbasename
            self._modalBase = self.mirrorModes
        elif mbasename == 'zonal':
            #print('Using zonal modes')
            self.modalBaseId = mbasename
            self._modalBase = self._createZonalMat()
        elif mbasename == 'hadamard':
            #print('Using Hadamard modes')
            self.modalBaseId = mbasename
            self._modalBase = self._createHadamardMat()
        else:
            #print('Using user-defined modes')
            self.modalBaseId = mbasename
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
        hadm = hadamard(2**numb)  # 892, 1 segment
        cmdBase = hadm[1:self._NActs+1, 1:self._NActs+1]
        #print('Removed 1st column of Hadamard matrix, or piston mode')
        return cmdBase
