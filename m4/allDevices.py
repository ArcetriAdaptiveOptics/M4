'''
@author: cs
'''

import numpy as np

#actCoord=

dimExt= 2.54
#dimInt=

        
class m4(object):
    def __init__(self):
        self._nActsTot= 5352
        
    def nActs(self):
        return self._nActsTot
    
    def pokeActs(self, indexing, amplitude, pushOrPull, modeMatrix= None):
        '''
        indexing= indicizzazione degli attuatori o dei modi
        pushOrPull= 1 o -1
        modeMatrix= nActsTot x nModes
        ''' 
        cmd= np.zeros(self._nActsTot)
        
        if modeMatrix is None:
            cmd[indexing]= amplitude * pushOrPull
            #comando per applicare cmd allo specchio 
        else:
            cmd=modeMatrix[:, indexing] * amplitude * pushOrPull
            #comando per applicare cmd allo specchio 
        


class segment(m4):
    def __init__(self, segmentIndex):
        super().__init__()
        
        if segmentIndex < 6:
            self._segmentIndex = segmentIndex
        else:
            raise OSError('Segment number %s doesnt exists' % segmentIndex)
        
        self._nActSeg= 892
        self._nSeg= 6
        
    def nActs(self):
        return self._nActSeg
    
    def pokeActs(self, indexing, amplitude, pushOrPull, modeMatrix= None):
        '''
        actuatorIndex= numero dell'attuatore del segmento (tra 0 e 799)
        '''
        if modeMatrix is None:
            act = indexing + (self._nActSeg * self._segmentIndex)
            super().pokeActs(act, amplitude, pushOrPull)
        else:
            super().pokeActs(indexing, amplitude, pushOrPull, modeMatrix)
            #così vorrò mettere 800 info in 5000 spazi quindi devo 
            #trovare come metterli nel posto giusto (con act) e tenere
            #gli altri a zero
        


