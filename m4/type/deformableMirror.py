'''
@author: cs
'''


class m4(object):
    def __init__(self):
        self._nActsTot= 5352
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
        
        self._nActSeg= 892
        self._nSeg= 6
        self._who= 'Segment number %s' % segmentIndex
        
    def nActs(self):
        return self._nActSeg