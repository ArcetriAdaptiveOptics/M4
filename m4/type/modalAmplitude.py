'''
@author: cs
'''


from m4.ground.configuration import Configuration
import os   
import pyfits
    
    
    
class ModalAmplitude():
    
    def __init__(self):
        self._modalAmplitude= None
        self._fitsfilename= None
    
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "ModalAmplitude")
    
    def getModalAmplitude(self):
        return self._modalAmplitude
    
    
    def save(self, tag, modalAmplitude):
        '''
            tag (stringa)= nome del file da salvare
            modalAmplitude (array)= vettore delle ampiezze
        '''
        storeInFolder= ModalAmplitude._storageFolder()
        filename= tag + '.fits'
        fitsFileName= os.path.join(storeInFolder, filename)
        pyfits.writeto(fitsFileName, modalAmplitude)
        
    
    @staticmethod 
    def load(fitsfilename):
        theObject= ModalAmplitude()
        storeInFolder= ModalAmplitude._storageFolder()
        allFitsFileName= os.path.join(storeInFolder, fitsfilename)
        hduList= pyfits.open(allFitsFileName)
        theObject._modalAmplitude= hduList[0].data
        theObject._fitsfilename= fitsfilename
        
        return theObject