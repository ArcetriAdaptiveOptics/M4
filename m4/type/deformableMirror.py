'''
@author: cs
'''

from m4.utils.configuration import Configuration

class m4(object):
    def __init__(self):
        self._nActsTot= Configuration.nActsTot
        self._who= 'All segments'
        
    def nActs(self):
        return self._nActsTot
    
    
class segment(m4):
    def __init__(self, segmentIndex):
        super().__init__()
        
        if segmentIndex < 6:
            self._segmentIndex = segmentIndex
        else:
            raise OSError('Segment number %s doesnt exists' % segmentIndex)
        
        self._nActSeg= Configuration.nActSeg
        self._nSeg= Configuration.n_Seg
        self._who= 'Segment number %s' % segmentIndex
        
    def nActs(self):
        return self._nActSeg