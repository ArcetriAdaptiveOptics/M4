'''
Authors
  - C. Selmi:  written in 2019
'''

import os
from m4.ground.timestamp import Timestamp

class TtFolder():
    """ Create a new folder using the generation of the tracking number"""

    def __init__(self, store_in_folder):
        """The constructor

        Parameters
        ----------
        store_in_folder: string
                        path where to create new fold
        """
        self._rootStoreFolder = store_in_folder


    def _createFolderToStoreMeasurements(self):
        """
        Returns
        -------
            dove: string
                all path generated
            tt: string
                only tracking number
        """
        self._tt = Timestamp.now()
        dove = os.path.join(self._rootStoreFolder, self._tt)
        if os.path.exists(dove):
            raise OSError('Directory %s exists' % dove)
        else:
            os.makedirs(dove)
        return dove, self._tt
