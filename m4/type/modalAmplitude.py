'''
@author: cs
'''


from m4.utils.configuration import Configuration
import os   
    
    
    
class ModalAmplitude():
    
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "ModalAmplitude")
    
    
    def save(self, tag):
        storeInFolder= ModalAmplitude._storageFolder()
        pass
    
    @staticmethod 
    def load(filename):
        pass