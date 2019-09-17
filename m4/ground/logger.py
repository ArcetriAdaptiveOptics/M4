'''
@author: cs
'''
import os

def log(info1, info2, info3, info4= None ):
    import logging
    
    if info4 is None:  
        string=" ".join([info1, info2, info3])
    else:
        string=" ".join([info1, info2, info3, info4])
    dove= '/Users/rm/Desktop/Arcetri/M4/ProvaCodice' 
    filename= os.path.join(dove, 'example.log')
         
    logging.basicConfig(filename= filename,
                        format='%(asctime)s - %(message)s')
    logging.warning(string)
        