'''
Authors
    - C. Selmi: written in 2021
'''

import os
import fnmatch
from m4.configuration.config import path_name

def ciccio(tt):
    rootPath = path_name.OPT_DATA_FOLDER
    pattern = tt

    for root, dirs, files in os.walk(rootPath):
        for directory in fnmatch.filter(dirs, pattern):
            final_path = os.path.join(root, directory)
    return final_path
