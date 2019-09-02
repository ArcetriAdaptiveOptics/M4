'''
@author: cs
'''

import os 
from m4.utils.configuration import Configuration
from m4.analyzerIFFunctions import AnalyzerIFF
from m4.utils.trackingNumberFolder import TtFolder
from m4.utils.tipTiltDetrend import TipTiltDetrend
import numpy as np
import pyfits


class Flattenig():
    
    def __init__(self, trackingNumber):
        self.iffPath= os.path.join(Configuration.FLATTENING_ROOT_FOLDER, trackingNumber)
        self._an= AnalyzerIFF.loadInfoFromh5Folder(self.iffPath)
        self._who= self._an._who
        self._command= None
        self._flatteningWf= None
        
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "Flattening")
        
        
    def flatCommand(self):
        #comando che permette di ottenere la misura del wf dall'interferometro (wf)
        wf=np.zeros(7)
        
        self._an.setDetectorMask(wf.mask | self._an.getMasterMask())
        rec= self._an.getReconstructor()
        wf_masked = np.ma.masked_array(wf.data, 
                                       mask=np.ma.mask_or(wf.mask, self._an.getMasterMask()))
        
        command= -np.dot(rec, wf_masked.compressed())
        
        return self._command
    
    def flattening(self, offset):
        #misuro la posizione dello specchio (pos)
        pos= np.zeros(7)
        
        cmd= pos + self._command
        #la applico e misuro il nuovo wf
        pass
    
    
    
    def save(self):
        storeInFolder= Flattenig._storageFolder()
        save= TtFolder(storeInFolder)
        dove, tt= save._createFolderToStoreMeasurements()
        fitsFileName= os.path.join(dove, 'info.fits')
        header= pyfits.Header()
        header['WHO']= self._who
        pyfits.writeto(fitsFileName, self._command, header)
        pyfits.append(fitsFileName, self._flatteningWf.data, header)
        pyfits.append(fitsFileName, self._flatteningWf.mask.astype(int), header)
        return tt
    
    
    @staticmethod
    def load(trackingNumber):
        theObject= Flattenig(trackingNumber)
        storeInFolder= Flattenig._storageFolder()
        folder= os.path.join(storeInFolder, trackingNumber)
        infoFitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.getheader(infoFitsFileName)
        hduList= pyfits.open(infoFitsFileName)
        theObject._who= header['WHO']
        theObject._command= hduList[0].data
        theObject._flattening= np.ma.masked_array(hduList[1].data, 
                                                  hduList[2].data.astype(bool))
        return theObject
            
            
            
            
            
        
        
        
        