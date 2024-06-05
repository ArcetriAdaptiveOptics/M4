#Preparazione - Class??
import os
import copy
import numpy as np
import h5py
from astropy.io import fits as pyfits
from m4.type.modalAmplitude import ModalAmplitude
from m4.type.modalBase import ModalBase
from m4.type.commandHisotry import CmdHistory
import configparser
config = configparser.ConfigParser() # Potrebbe ancora servire, vediamo

class IFFCapturePreparation():

    def __init__(self): #file .ini per inizializzare? No, non c'è bisogno. leggono le funzioni in caso
        self._nActs         = # get N Acts from file? Or load DM?
        self._cmdMatrix     = # load mirror command matrix, Serve?
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
        # Gli input sono gli output di altre due funzioni (interne?)
        # Ha senso? Ni, dipende da quanto sono complesse queste due funzioni
        # e quanto separarle rende il codice flessibile e chiaro. Sicuramente 
        # il salvataggio o il loading va incluso

        return aux_cmdHistory

    # Ha senso tenere le due funzioni separate?

    def createRegistrationPattern(self, reg_modes, reg_amplitudes, reg_template):
        # Sono dei modi sparati in successione, tipo due mega trifogli, che dicono
        # 'da qui iniziano i dati'. Sono gli Rx e Ry.

        return reg_patternCmdList

    def createTriggerPadding(self, padding_template=None, trigger_mode=None, trigger_amplitude=None, trigger_template=None, save=False):# Salvare gli output? Vengono salvati comunque nella classe

        # Come vogliamo gli imput? File? Vettori passati da python? Info su file di configurazione? 
        # Catena di if per le variabili di ingresso. Se l'imput è None, allora carica le variabili
        # 'standard', altrimenti vengono specificate

        # Serve ad identificare il primo frame buono per far partire l'acquisizione / unwrapping. E' il T_0

        # Ok no va bene perché c'è la possibilità che vogliano essere cambiati o che si vogliano skippare
        # quindi in generale è meglio fare le funzioni. Gli imput vengono da file. Caricato all'inizio? forse meglio.

        return trigger_paddingCmdList


