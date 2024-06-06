#Preparazione - Class??
import os
import copy
import numpy as np
import h5py
from astropy.io import fits as pyfits
import configparser
from m4.mini_ott import timehisory as th
from m4.configuration import iffConfig as iffc  #contains a dictionary or similar with configuration

config = configparser.ConfigParser() # Potrebbe ancora servire, vediamo

class IFFCapturePreparation():

    def __init__(self): #file .ini per inizializzare? No, non c'è bisogno. leggono le funzioni in caso
        self._nActs         = # get N Acts from file? Or load DM?
        self._cmdMatrix     = # loaded command base in _getCmdMatrix()
        self._indexingList  = None

        self._paddingTemp   = None
        self._triggMode     = None
        self._triggAmpl     = None
        self._triggTemp     = None
        self._regModes      = None
        self._regAmpl       = None
        self._regTemp       = None
        
        self.triggPadCmdList   = None
        self.regPadCmdList     = None
        self.auxCmdHistory     = None
        self.cmdMatHistory     = None


    def createTimedCmdHistory(self, aux_cmdhistory, cmd_matrixHistory, timing_info): 
        # timing info as dict? ini?? -> comunque esterne al codice
        # Si può prendere spunto dalle funzioni sotto

        return timedCmdList

    def createCmdMatrixHistory(self, modes_matrix, modes_amplitude, cmd_template, shuffle=False):
        # Si può riutilizzare molto di quello che si ha già, 
        # magari sistemato un po' e reso più flessibile e chiaro
        self._modeVector = copy.copy(mode_vector) # modes_matrix?
        self._nPushPull = n_push_pull # ??
        self._cmdMatrix = cmd_matrix # cmd_template?
        indList = []

        if shuffle==True:
            for i in range(len(n_push_pull):
                np.random.shuffle(mode_vector)
                indexingList.append(list(mode_vector))
        else:
            for i in range(n_push_pull):
                indList.append(mode_vector)

        self._indexingList = np.array(indList)

        n_frame = mode_vector.size * n_push_pull
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

    def createAuxCmdHistory(self, reg_patterncmdList, trigger_paddingCmdList):
    '''
    Creates the initial parteof the final command history matrix that will be passed to M4. This includes the Trigger Frame, the first frame to have a non-zero command, and the Padding Frame, two frames with high rms, useful for setting a start to the real acquisition.

    Parameters
    ----------


    Result
    ------
    '''
        # Gli input sono gli output di altre due funzioni (interne?)
        # Ha senso? Ni, dipende da quanto sono complesse queste due funzioni
        # e quanto separarle rende il codice flessibile e chiaro. Sicuramente 
        # il salvataggio o il loading va incluso
        mette insieme i due sotto
        return aux_cmdHistory

    # Ha senso tenere le due funzioni sotto separate da questa sopra? Potrebbee non essere necessario

    def createRegistrationPattern(self, reg_modes, reg_amplitudes, reg_template):
        # Sono dei modi sparati in successione, tipo due mega trifogli, che dicono
        # 'da qui iniziano i dati'. Sono gli Rx e Ry.
        stesso schema di trigger padding
        input viene da iffc: e sarebbe quali attuatori alzare per fare i marker, con che ampiezza e eventualmente quante ripetizioni (template)
        return reg_patternCmdList

    def createTriggerPadding(self, padding_template=None, trigger_mode=None, trigger_amplitude=None, trigger_template=None, save=False):# Salvare gli output? Vengono salvati comunque nella classe

        # Come vogliamo gli imput? File? Vettori passati da python? Info su file di configurazione? 
        # Catena di if per le variabili di ingresso. Se l'imput è None, allora carica le variabili
        # 'standard', altrimenti vengono specificate

        # Serve ad identificare il primo frame buono per far partire l'acquisizione / unwrapping. E' il T_0

        # Ok no va bene perché c'è la possibilità che vogliano essere cambiati o che si vogliano skippare
        # quindi in generale è meglio fare le funzioni. Gli imput vengono da file. Caricato all'inizio? forse meglio.

        cmdxx= getCmdMatrix(mirrorModes) #legge da configurazioni
        triggCmd = cmdxx[triggId,:]
        paddingscheme = lo leggi da iffc (vedi import iniziale, sono quanti zeri mettere all'inizio)
        triggTimeHist = la creai asemblando gli zeri con i lcomando di trigg



        return trigger_paddingCmdList


    def _getCmdMatrix(identif, mlist):
        """
        This function ...

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
