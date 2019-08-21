'''
@author: cs
'''

import os
import pyfits
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
    
    
    def _saveIFFInfoModal(self, folder, who, amplitude,
                                vectorOfActuatorsOrModes, nPushPull,
                                 indexingList, cmdMatrix):
            
        fitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.Header()
        header['NPUSHPUL']= nPushPull
        header['WHO']= who
        pyfits.writeto(fitsFileName, vectorOfActuatorsOrModes, header)
        pyfits.append(fitsFileName, indexingList, header)
        pyfits.append(fitsFileName, cmdMatrix, header)
        pyfits.append(fitsFileName, amplitude, header)
        
        
    def _saveIFFInfoZonal(self, folder, who, amplitude,
                                vectorOfActuatorsOrModes, nPushPull):
            
        fitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.Header()
        header['NPUSHPUL']= nPushPull
        header['WHO']= who
        pyfits.writeto(fitsFileName, vectorOfActuatorsOrModes, header)
        pyfits.append(fitsFileName, amplitude, header)
            
            
    @staticmethod
    def loadIFFInfoModal(folder):
        additionalInfoFitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.getheader(additionalInfoFitsFileName)
        hduList= pyfits.open(additionalInfoFitsFileName)
        actsVector= hduList[0].data
        indexingList= hduList[1].data
        cmdMatrix= hduList[2].data
        cmdAmpl= hduList[3].data
        who= header['WHO']
        try:
            nPushPull= header['NPUSHPUL']
        except KeyError:
            nPushPull= 1
            
        return (who, actsVector, cmdMatrix, indexingList, cmdAmpl, nPushPull)
    
    @staticmethod
    def loadIFFInfoZonal(folder):
        additionalInfoFitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.getheader(additionalInfoFitsFileName)
        hduList= pyfits.open(additionalInfoFitsFileName)
        actsVector= hduList[0].data
        cmdAmpl= hduList[1].data
        who= header['WHO']
        try:
            nPushPull= header['NPUSHPUL']
        except KeyError:
            nPushPull= 1
            
        return (who, actsVector, cmdAmpl, nPushPull)
    