#Preparazione - Class??
from astropy.io import fits as pyfits
from m4.type.modalAmplitude import ModalAmplitude
from m4.type.modalBase import ModalBase
from m4.type.commandHisotry import CmdHistory
import configparser
config = configparser.ConfigParser() # Potrebbe ancora servire, vediamo

class Preparation():

    def __init__(self): #file .ini per inizializzare? No, non c'è bisogno. leggono le funzioni in caso
        self._nActs         = # get N Acts from file? Or load DM?
        self._cmdMatrix     = # load mirror command matrix, Serve?

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

    def createCmdMatrix(self, modes_matrix, modes_amplitude, cmd_template, shuffle=False):
        # Si può riutilizzare molto di quello che si ha già, 
        # magari sistemato un po' e reso più flessibile e chiaro
        

        return cmd_matrixHistory

    def createAuxCmdHistory(self, reg_patterncmdList, trigger_paddingCmdList):
        # Gli input sono gli output di altre due funzioni (interne?)
        # Ha senso? Ni, dipende da quanto sono complesse queste due funzioni
        # e quanto separarle rende il codice flessibile e chiaro. Sicuramente 
        # il salvataggio o il loading va incluso

        return aux_cmdHistory

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


