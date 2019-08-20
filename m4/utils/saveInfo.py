'''
@author: cs
'''

import os
import pyfits
import numpy as np
from m4.utils.timestamp import Timestamp

class SaveAdditionalInfo(object):
    
    def __init__(self, storeInFolder):
        self._rootStoreFolder= storeInFolder
        
    
    def _createFolderToStoreMeasurements(self):
        self._tt=Timestamp.now()
        dove= os.path.join(self._rootStoreFolder, self._tt)
        if os.path.exists(dove):
            raise OSError('Directory %s exists' % dove)
        else:
            os.makedirs(dove)
        return dove
    
    
    def _saveIFFInfo(self, folder, who, amplitude,
                                vectorOfActuatorsOrModes, nPushPull,
                                 indexingList= None, cmdMatrix= None):
#         if (indexingList, cmdMatrix) is None:
#             indexingList=np.zeros(3)
#             cmdMatrix= np.zeros((3,3))
#         else:
#             indexingList=indexingList
#             cmdMatrix=cmdMatrix
            
        fitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.Header()
        header['CMDAMPL']= amplitude
        header['NPUSHPUL']= nPushPull
        header['WHO']= who
        pyfits.writeto(fitsFileName, vectorOfActuatorsOrModes, header)
        pyfits.append(fitsFileName, indexingList, header)
        pyfits.append(fitsFileName, cmdMatrix, header)
            
            
    @staticmethod
    def loadAdditionalInfo(folder):
        additionalInfoFitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.getheader(additionalInfoFitsFileName)
        hduList= pyfits.open(additionalInfoFitsFileName)
        actsVector= hduList[0].data
        indexingList= hduList[1].data
        cmdMatrix= hduList[2].data
        cmdAmpl= header['CMDAMPL']
        who= header['WHO']
        try:
            nPushPull= header['NPUSHPUL']
        except KeyError:
            nPushPull= 1
            
        return (who, actsVector, cmdMatrix, indexingList, cmdAmpl, nPushPull)