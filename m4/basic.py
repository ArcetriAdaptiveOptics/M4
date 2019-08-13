'''
@author: cs
'''

import m4 as _m4
import os
import pyfits


class SaveInfo(object):
    
    def log(self, info1, info2):
        import logging
        
        string=" ".join([info1, info2])
         
        logging.basicConfig(filename='example.log',format='%(asctime)s - %(message)s')
        logging.warning(string)
        
    def _saveIFFInfo(self, folder, amplitude,
                            vectorOfActuatorsOrModes, nPushPull):
        fitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.Header()
        header['CMDAMPL']= amplitude
        header['NPUSHPUL']= nPushPull
        pyfits.writeto(fitsFileName, vectorOfActuatorsOrModes, header)


class PhaseMap(object):
    
    def absPhaseMap(self):
        pass

    def diffPhaseMap(self):
        pass
    
    pass



class ROI(object):
    
    def selectMasterRoi(self):
        pass
    
    def selectSecondaryRoi(self):
        pass
    
    pass


class ZernikeGenerator(object):
    
    pass




    
    