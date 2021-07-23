'''
Authors
  - C. Selmi: written in 2019
             modified in July 2021
'''

import os
import fnmatch
from m4.ground.timestamp import Timestamp
from m4.configuration.config import path_name


def createFolderToStoreMeasurements(store_in_folder):
    """
    Create a new folder using the generation of the tracking number
    Parameters
    ----------
    store_in_folder: string
                    path where to create new fold
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

def findTrackingNumberPath(tt):
    ''' Function to find the entire path from tracking number
    Parameters
    ----------
        tt: string
            tracking number
    Returns
    -------
        final_path: string
            total path
    '''
    rootPath = path_name.OPT_DATA_FOLDER
    pattern = tt

    for root, dirs, files in os.walk(rootPath):
        for directory in fnmatch.filter(dirs, pattern):
            final_path = os.path.join(root, directory)
    try:
        return final_path
    except:
        raise OSError('Tracking number %s does not exists' % tt)

# '''
# Authors
#   - C. Selmi:  written in 2019
# '''
# 
# import os
# from m4.ground.timestamp import Timestamp
# 
# 
# 
# def makeTrackingNumberFolder(store_in_folder):
#     """
#     Parameters
#     ----------
#     store_in_folder:string
#         path location for new tracking number folder
# 
#     Returns
#     -------
#         dove: string
#             all path generated
#         tt: string
#             only tracking number
#     """
#     tt = Timestamp.now()
#     dove = os.path.join(store_in_folder, tt)
#     if os.path.exists(dove):
#         raise OSError('Directory %s exists' % dove)
#     else:
#         os.makedirs(dove)
#     return dove, tt
