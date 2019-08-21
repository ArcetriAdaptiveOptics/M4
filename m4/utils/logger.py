'''
@author: cs
'''
import os

def log(info1, info2):
    import logging
        
    string=" ".join([info1, info2])
    dove= '/Users/rm/Desktop/Arcetri/M4/ProvaCodice' 
    filename= os.path.join(dove, 'example.log')
         
    logging.basicConfig(filename= filename,
                        format='%(asctime)s - %(message)s')
    logging.warning(string)
        