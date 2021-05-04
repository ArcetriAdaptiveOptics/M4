'''
Authors
  - C. Selmi:  written in 2019
'''

import os
from m4.ground.timestamp import Timestamp



def makeTrackingNumberFolder(store_in_folder):
    """
    Parameters
    ----------
    store_in_folder:string
        path location for new tracking number folder

    Returns
    -------
        dove: string
            all path generated
        tt: string
            only tracking number
    """
    tt = Timestamp.now()
    dove = os.path.join(store_in_folder, tt)
    if os.path.exists(dove):
        raise OSError('Directory %s exists' % dove)
    else:
        os.makedirs(dove)
    return dove, tt
