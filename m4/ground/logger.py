'''
@author: cs
'''
import os
from m4.ground.configuration import Configuration

def log(info1, info2, info3, info4= None ):
    import logging
    
    if info4 is None:  
        string=" ".join([info1, info2, info3])
    else:
        string=" ".join([info1, info2, info3, info4])
    dove= Configuration.LOG_ROOT_FOLDER
    filename= os.path.join(dove, 'example.log')
         
    logging.basicConfig(filename= filename,
                        format='%(asctime)s - %(message)s')
    logging.warning(string)
        