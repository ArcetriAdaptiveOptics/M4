'''
Authors
  - C. Selmi: written in 2022
  - L. Oggioni: added some functionalities
'''

import os
import logging
import numpy as np
from m4.devices.base_deformable_mirror import BaseDeformableMirror
from m4.configuration.ott_parameters import M4Parameters
from m4.configuration import config_folder_names as conf
from m4.ground import read_data


class FakeM4DM(BaseDeformableMirror):
    '''
    '''

    def __init__(self):
        """The constructor """
        self._actPos = np.zeros(M4Parameters.N_ACT_SEG)
        self._logger = logging.getLogger('FakeM4DM')
        self.m4pupil = read_data.readFits_data(os.path.join(conf.MIRROR_FOLDER,
                                                        conf.mirror_conf,
                                                        'm4_mech_pupil-bin2.fits'))
        self.m4ima = self.m4pupil * 0.
        self.CapsensGain = np.load(os.path.join(conf.MIRROR_FOLDER,
                                                conf.mirror_conf, 'CapsensGain.npy'))
        self.ifmat = read_data.readFits_data(os.path.join(conf.MIRROR_FOLDER,
                                                    conf.mirror_conf,
                                                   'if_sect4_rot-bin2.fits'))
        self.ifidx = read_data.readFits_data(os.path.join(conf.MIRROR_FOLDER,
                                                           conf.mirror_conf,
                                                    'if_idx4_rot-bin2.fits'))
        self.vmat = read_data.readFits_data(os.path.join(conf.MIRROR_FOLDER,
                                                   conf.mirror_conf, 'ff_v_matrix.fits'))

    def setActsCommand(self, command, delta=True):
        self._actPos = command
        image = self._mirrorCommand(self._actPos, delta)
        return

    def getActsCommand(self):
        return self._actPos

    def getNActs(self):
        return M4Parameters.N_ACT_SEG
    
    
    def act_incr(self, inc, rel=False):
    
        ''' comando di posizione a tutti i pistoni 
        
        Parameters
        ----------
        inc: float
            incremento in m
        rel: logic
            incremento relativo (True) o posizione assoluta (False)
        
        '''
        
        inc=inc*1e6
        comm=np.ones(892)*inc
        self.setActsCommand(comm,rel)
        return
        
    def act_zero(self):
        
        ''' azzera la posizione dei pistoni   
        '''      
        comm=np.zeros(892)
        rel=False
        self.setActsCommand(comm,rel)
        return
        
    def act_random(self, ampiezza,rel=False):
        
        ''' genera una distribuzione casuale di pistonate tra 0 e "ampiezza"  
                 
        Parameters
        ----------
        
        ampiezza: float
            massima ampiezza (m) generata dalla distribuzione casuale di pistoni 
        
        rel: logic
            incremento relativo (True) o posizione assoluta (False)
        
        '''    
        ampiezza=ampiezza*1e6
        comm=np.random.rand(892)*ampiezza
        self.setActsCommand(comm,rel)
        return

    def _generateAndSaveCapsensGain(self, file_name):
        ''' Function to generate and save a random gain 
        Parameters
        ----------
        file_name: string
            cap sens gain file name (.npy)

        Returns
        ------
        file_path: string
            cap sens complete file path
        '''
        gain = np.ones(self.getNActs()) + np.random.rand(self.getNActs())*0.1 - 0.05
        file_path = os.path.join(conf.MIRROR_FOLDER, conf.mirror_conf, file_name)
        np.save(file_path, gain)
        return file_path

    def _applyCapsensGain(self, command):
        '''
        Parameters
        ----------
        command: numpy array [892]
            vector for multiply capsens gain (for real effect)

        Returns
        ------
        comm: numpy array
            vector multiplied
        '''

        command = command*self.CapsensGain
        return command

    def _mirrorCommand(self, comm, delta=True):
        '''
        Parameters
        ----------
        comm: numpy array [892]
            command for a segment

        Other Parameters
        ----------------
        delta: boolean
            if delta is True relative command is used
            else absolute command

        Return
        ------
        self.m4ima: numpy array
            m4 total image
        '''

        comm = self._applyCapsensGain(comm)

        if delta==True:
            self.m4ima.flat[self.ifidx] = self.m4ima.flat[self.ifidx]+np.matmul(np.transpose(self.ifmat), comm)
        elif delta==False:
            self.m4ima.flat[self.ifidx] = np.matmul(np.transpose(self.ifmat), comm)
        return self.m4ima

