'''
@author: cs
'''

import os
from m4.utils.timestamp import Timestamp

class TtFolder(object):
    
    def __init__(self, storeInFolder):
        self._rootStoreFolder= storeInFolder
        
    
    def _createFolderToStoreMeasurements(self):
        self._tt=Timestamp.now()
        dove= os.path.join(self._rootStoreFolder, self._tt)
        if os.path.exists(dove):
            raise OSError('Directory %s exists' % dove)
        else:
            os.makedirs(dove)
        return dove, self._tt

    

    